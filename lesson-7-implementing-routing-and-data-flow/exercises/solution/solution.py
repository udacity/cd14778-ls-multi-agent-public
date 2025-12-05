import os
import time
from collections import deque
from typing import Deque, Dict, List, Optional, Union

import dotenv
from langchain.agents import AgentExecutor, create_openai_tools_agent
from langchain_core.prompts import ChatPromptTemplate, MessagesPlaceholder
from langchain_core.tools import tool
from langchain_openai import ChatOpenAI

# --- Configuration and Setup ---

dotenv.load_dotenv("../../.env")
OPENAI_API_KEY = os.getenv("OPENAI_API_KEY")

# NOTE: Bluesky is optional. If you'd like to try on real data from the
#       Bluesky firehose, set to True. Otherwise simulated data will be used
USE_BLUESKY_DATA = False
BLUESKY_POST_LIMIT = 10 # post limit to ingest (if using Bluesky) - keep below 100 for performance

try:
    from atproto import FirehoseSubscribeReposClient, firehose_models, models
except Exception:
    print("atproto not installed; falling back to synthetic data")

ENGAGEMENT_THRESHOLD = 1000

llm = ChatOpenAI(
    model="gpt-4o-mini",
    api_key=OPENAI_API_KEY,
    base_url="https://openai.vocareum.com/v1",
    temperature=0,
)

Post = Dict[str, Union[str, int]]

# --- State Management: The Output Queues ---

class FactCheckQueues:
    def __init__(self):
        self.vaccine_safety_queue: Deque[Post] = deque()
        self.treatment_efficacy_queue: Deque[Post] = deque()

    def add_to_vaccine_queue(self, post: Post, is_priority: bool = False):
        if is_priority:
            self.vaccine_safety_queue.appendleft(post)
            print("  🔥 PRIORITY: Post added to the FRONT of the vaccine queue.")
        else:
            self.vaccine_safety_queue.append(post)
            print("  ✅ Post added to the vaccine queue.")

    def add_to_treatment_queue(self, post: Post, is_priority: bool = False):
        if is_priority:
            self.treatment_efficacy_queue.appendleft(post)
            print("  🔥 PRIORITY: Post added to the FRONT of the treatment queue.")
        else:
            self.treatment_efficacy_queue.append(post)
            print("  ✅ Post added to the treatment queue.")

    def display_queues(self):
        print("\n" + "="*80)
        print("FINAL QUEUE STATUS")
        print("="*80)
        print("\n🔬 Vaccine Safety Fact-Check Queue:")
        if not self.vaccine_safety_queue:
            print("  [Empty]")
        else:
            for i, post in enumerate(self.vaccine_safety_queue):
                priority_tag = "🔥" if post["engagement"] > ENGAGEMENT_THRESHOLD else " "
                print(f"  {i+1}. {priority_tag} Post ID: {post['id']} (Engagement: {post['engagement']})")

        print("\n💊 Treatment Efficacy Fact-Check Queue:")
        if not self.treatment_efficacy_queue:
            print("  [Empty]")
        else:
            for i, post in enumerate(self.treatment_efficacy_queue):
                priority_tag = "🔥" if post["engagement"] > ENGAGEMENT_THRESHOLD else " "
                print(f"  {i+1}. {priority_tag} Post ID: {post['id']} (Engagement: {post['engagement']})")
        print("\n" + "="*80)
        print("✅ SUCCESS: The system demonstrated content-based classification and priority-based routing.")
        print("="*80)

fact_check_queues = FactCheckQueues()

# --- Tools for the Router Agent (STRICT: require both args; no defaults) ---

@tool("route_to_vaccine_fact_check")
def route_to_vaccine_fact_check(post: Post, is_priority: bool) -> str:
    """
    Route a vaccine-safety claim.
    Args:
      post: dict with keys {'id': str, 'text': str, 'engagement': int}
      is_priority: boolean
    """
    if not isinstance(post, dict) or not {"id", "text", "engagement"} <= set(post.keys()):
        raise ValueError("Bad 'post' payload; expected dict with keys id, text, engagement.")
    fact_check_queues.add_to_vaccine_queue(post, is_priority)
    return "Post successfully routed to vaccine safety queue."

@tool("route_to_treatment_fact_check")
def route_to_treatment_fact_check(post: Post, is_priority: bool) -> str:
    """
    Route a treatment-efficacy claim.
    Args:
      post: dict with keys {'id': str, 'text': str, 'engagement': int}
      is_priority: boolean
    """
    if not isinstance(post, dict) or not {"id", "text", "engagement"} <= set(post.keys()):
        raise ValueError("Bad 'post' payload; expected dict with keys id, text, engagement.")
    fact_check_queues.add_to_treatment_queue(post, is_priority)
    return "Post successfully routed to treatment efficacy queue."

# --- Data Source ---

def get_social_media_firehose(use_simulated: bool = True, post_limit: int = BLUESKY_POST_LIMIT) -> List[Post]:
    if use_simulated:
        print("🔥 Using SIMULATED social media firehose...")
        return [
            {'id': 'post_001', 'text': 'Heard the new X-4 antibody treatment completely neutralizes the virus in hours.', 'engagement': 250},
            {'id': 'post_002', 'text': 'My arm was so sore after the Vax-Protect shot. I think it might be causing long-term nerve damage.', 'engagement': 800},
            {'id': 'post_003', 'text': 'This is huge! Dr. Sharma says the Z-MAB drug has a 95% success rate. Retweet this!', 'engagement': 15000},
            {'id': 'post_004', 'text': 'Just saw a report that the government is hiding data about adverse reactions to the mRNA vaccine.', 'engagement': 5200},
            {'id': 'post_005', 'text': 'Why is no one talking about herbal remedies? A simple tea can stop the pathogen!', 'engagement': 150},
            {'id': 'post_006', 'text': 'If you feel sick after your shot, it just means it\'s working!', 'engagement': 30},
            {'id': 'post_007', 'text': 'Source says a new pill, "CureAll", is being blocked by big pharma but it\'s the real deal.', 'engagement': 950},
        ]
    else:
        print(f"🔌 Connecting to LIVE Bluesky firehose. Waiting for {post_limit} new posts...")

        client = FirehoseSubscribeReposClient()
        posts: List[Post] = []

        def on_message_handler(message: 'firehose_models.MessageFrame') -> None:
            if not isinstance(message.data, firehose_models.Commit):
                return
            if not message.data.blocks or not message.data.ops:
                return

            for op in message.data.ops:
                if op.action == 'create' and op.path.endswith('app.bsky.feed.post'):
                    record = models.get_model_as_dict(
                        models.load_data(message.data.blocks, op.cid)
                    )

                    if record.get('text'):
                        formatted_post: Post = {
                            # 🔧 safer CID stringification
                            'id': str(op.cid),
                            'text': record.get('text', ''),
                            'engagement': 0
                        }
                        posts.append(formatted_post)
                        print(f"   ({len(posts)}/{post_limit}) Post captured from firehose...")
                        if len(posts) >= post_limit:
                            client.stop()

        client.start(on_message_handler)
        return posts

# --- Main Orchestration Logic ---

def main(use_simulated_firehose: bool):
    # 1) Classifier (simple chain)
    classifier_prompt = ChatPromptTemplate.from_messages([
        ("system",
         "You are an expert content analyst. Classify the social media post into "
         "'vaccine_safety_claim' or 'treatment_efficacy_claim'. "
         "Respond with ONLY the category name."),
        ("user", "{post_text}"),
    ])
    classifier_chain = classifier_prompt | llm

    # 2) Router agent (static tools, strict, with self-correction)
    router_tools = [route_to_vaccine_fact_check, route_to_treatment_fact_check]

    # Surface tool errors back to the LLM as observations so it can retry
    for t in router_tools:
        t.handle_tool_error = True

    router_prompt = ChatPromptTemplate.from_messages([
        ("system",
         "You are an intelligent router. Choose exactly ONE tool.\n"
         "Routing policy:\n"
         "- If claim_type == 'vaccine_safety_claim' -> call route_to_vaccine_fact_check.\n"
         "- If claim_type == 'treatment_efficacy_claim' -> call route_to_treatment_fact_check.\n"
         "You MUST pass BOTH arguments to the selected tool: "
         "post (dict exactly as provided) and is_priority (boolean). "
         "Do not rename keys or change values."),
        ("human",
         "Routing inputs:\n"
         "claim_type: {claim_type}\n"
         "is_priority: {is_priority}\n"
         "post: {post}"),
        MessagesPlaceholder("agent_scratchpad"),
    ])

    router_agent = create_openai_tools_agent(llm, router_tools, router_prompt)

    router_executor = AgentExecutor(
        agent=router_agent,
        tools=router_tools,
        verbose=True,
        handle_parsing_errors=True,  # <- let the agent recover if it formats wrong / missing keys
        max_iterations=3,             # <- avoid infinite retries
    )

    # 3) Get data and process
    posts = get_social_media_firehose(use_simulated=use_simulated_firehose)

    global ENGAGEMENT_THRESHOLD
    if not use_simulated_firehose:
        ENGAGEMENT_THRESHOLD = 0
        print(f"\n⚠️ Live mode active. Engagement threshold lowered to {ENGAGEMENT_THRESHOLD}.")

    for i, post in enumerate(posts):
        print("\n" + "-"*60)
        print(f"Processing Post {i+1}/{len(posts)} (ID: {post['id']})")
        print(f'Text: "{post["text"]}"')
        print(f"Engagement: {post['engagement']}")

        # Priority
        is_priority = post["engagement"] > ENGAGEMENT_THRESHOLD
        print(f"  🚦 Priority Check: {'HIGH' if is_priority else 'Normal'} (>{ENGAGEMENT_THRESHOLD} engagement)")

        # Classify
        print("  🧠 Classifying content with Classifier Agent...")
        classification_result = classifier_chain.invoke({"post_text": post["text"]})
        claim_type_raw = classification_result.content.strip().lower()
        claim_type = "vaccine_safety_claim" if "vaccine" in claim_type_raw else "treatment_efficacy_claim"
        print(f"  - Agent Classification: {claim_type} (raw: {claim_type_raw})")

        # Route (agentic tool call with self-correction)
        print("  - Routing to queue with Router Agent...")
        router_executor.invoke({
            "input": "Route now.",
            "post": post,                 # required; if LLM omits it, it will see the error and retry
            "is_priority": is_priority,   # required
            "claim_type": claim_type,     # guides tool choice
        })

        time.sleep(0.5)

    # 4) Show results
    fact_check_queues.display_queues()

if __name__ == "__main__":
    main(use_simulated_firehose=not USE_BLUESKY_DATA)
