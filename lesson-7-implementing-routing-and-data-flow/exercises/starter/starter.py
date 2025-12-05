"""
-- INSTRUCTIONS (brief) --
Finish the three TODO sections in main():
  1) Build the classifier chain (prompt -> LLM) that outputs ONLY:
     'vaccine_safety_claim' or 'treatment_efficacy_claim'.
  2) Compute priority and get the claim_type from the classifier.
  3) Use the router agent to call exactly ONE tool with BOTH args:
     post (dict) and is_priority (bool).

Bluesky live data is optional. Default uses simulated posts.
"""

import os
import time
from collections import deque
from typing import Deque, Dict, List, Union

import dotenv
from langchain.agents import AgentExecutor, create_openai_tools_agent
from langchain_core.prompts import ChatPromptTemplate, MessagesPlaceholder
from langchain_core.tools import tool
from langchain_openai import ChatOpenAI

# ========== Configuration ==========
dotenv.load_dotenv("../../.env")
OPENAI_API_KEY = os.getenv("OPENAI_API_KEY")

USE_BLUESKY_DATA = False
BLUESKY_POST_LIMIT = 10

# Try to enable Bluesky only if requested
try:
    if USE_BLUESKY_DATA:
        from atproto import FirehoseSubscribeReposClient, firehose_models, models
except Exception:
    print("ℹ️  'atproto' not available; falling back to simulated data.")
    USE_BLUESKY_DATA = False

ENGAGEMENT_THRESHOLD = 1000 if not USE_BLUESKY_DATA else 0

llm = ChatOpenAI(
    model="gpt-4o-mini",
    api_key=OPENAI_API_KEY,
    base_url="https://openai.vocareum.com/v1",
    temperature=0,
)

Post = Dict[str, Union[str, int]]

# ========== State (Queues) ==========
class FactCheckQueues:
    def __init__(self):
        self.vaccine_safety_queue: Deque[Post] = deque()
        self.treatment_efficacy_queue: Deque[Post] = deque()

    def add_to_vaccine_queue(self, post: Post, is_priority: bool = False):
        print(f"  ➡️  Route → VACCINE (Priority: {is_priority})")
        self.vaccine_safety_queue.appendleft(post) if is_priority else self.vaccine_safety_queue.append(post)

    def add_to_treatment_queue(self, post: Post, is_priority: bool = False):
        print(f"  ➡️  Route → TREATMENT (Priority: {is_priority})")
        self.treatment_efficacy_queue.appendleft(post) if is_priority else self.treatment_efficacy_queue.append(post)

    def display_queues(self):
        print("\n" + "="*80 + "\n✅ FINAL QUEUE STATUS\n" + "="*80)
        print("\n🔬 Vaccine Safety:")
        if not self.vaccine_safety_queue: print("  [Empty]")
        for i, p in enumerate(self.vaccine_safety_queue):
            tag = "🔥" if p["engagement"] > ENGAGEMENT_THRESHOLD else " "
            print(f"  {i+1}. {tag} ID: {p['id']} (eng: {p['engagement']})")
        print("\n💊 Treatment Efficacy:")
        if not self.treatment_efficacy_queue: print("  [Empty]")
        for i, p in enumerate(self.treatment_efficacy_queue):
            tag = "🔥" if p["engagement"] > ENGAGEMENT_THRESHOLD else " "
            print(f"  {i+1}. {tag} ID: {p['id']} (eng: {p['engagement']})")
        print("\n" + "="*80)

fact_check_queues = FactCheckQueues()

# ========== Tools ==========
@tool("route_to_vaccine_fact_check")
def route_to_vaccine_fact_check(post: Post, is_priority: bool):
    """Route a vaccine-safety claim. Args: post (dict), is_priority (bool)."""
    if not isinstance(post, dict) or not {"id","text","engagement"} <= set(post.keys()):
        raise ValueError("Bad 'post' payload; expected dict with keys id, text, engagement.")
    fact_check_queues.add_to_vaccine_queue(post, is_priority)
    return "Post routed to vaccine safety."

@tool("route_to_treatment_fact_check")
def route_to_treatment_fact_check(post: Post, is_priority: bool):
    """Route a treatment-efficacy claim. Args: post (dict), is_priority (bool)."""
    if not isinstance(post, dict) or not {"id","text","engagement"} <= set(post.keys()):
        raise ValueError("Bad 'post' payload; expected dict with keys id, text, engagement.")
    fact_check_queues.add_to_treatment_queue(post, is_priority)
    return "Post routed to treatment efficacy."

# ========== Data ==========
def get_simulated_firehose() -> List[Post]:
    print("🔥 Using SIMULATED social media firehose...")
    return [
        {'id':'post_001','text':'Heard the new X-4 antibody treatment completely neutralizes the virus in hours.','engagement':250},
        {'id':'post_002','text':'My arm was so sore after the Vax-Protect shot. I think it might be causing long-term nerve damage.','engagement':800},
        {'id':'post_003','text':'This is huge! Dr. Sharma says the Z-MAB drug has a 95% success rate. Retweet this!','engagement':15000},
        {'id':'post_004','text':'Just saw a report that the government is hiding data about adverse reactions to the mRNA vaccine.','engagement':5200},
        {'id':'post_005','text':'Why is no one talking about herbal remedies? A simple tea can stop the pathogen!','engagement':150},
    ]

def get_bluesky_firehose(limit: int = BLUESKY_POST_LIMIT) -> List[Post]:
    print(f"🔌 Connecting to LIVE Bluesky firehose… collecting up to {limit} posts.")
    client = FirehoseSubscribeReposClient()
    posts: List[Post] = []

    def on_message_handler(message: 'firehose_models.MessageFrame') -> None:
        if not isinstance(message.data, firehose_models.Commit):
            return
        if not message.data.blocks or not message.data.ops:
            return
        for op in message.data.ops:
            if op.action == "create" and op.path.endswith("app.bsky.feed.post"):
                record = models.get_model_as_dict(models.load_data(message.data.blocks, op.cid))
                text = record.get("text", "")
                if text:
                    posts.append({'id': str(op.cid), 'text': text, 'engagement': 0})
                    print(f"  ({len(posts)}/{limit}) captured…")
                    if len(posts) >= limit:
                        client.stop()

    client.start(on_message_handler)
    return posts

def get_firehose(use_bluesky: bool) -> List[Post]:
    return get_bluesky_firehose() if use_bluesky else get_simulated_firehose()

# ========== Main (Student Work) ==========
def main():
    posts = get_firehose(USE_BLUESKY_DATA)

    # -------- STUDENT TODO 1: DEFINE YOUR AGENTS --------
    # (a) Classifier: prompt -> llm chain that returns ONLY the label.
    classifier_prompt = ChatPromptTemplate.from_messages([
        # system message goes here
        # must instruct the model to output ONLY:
        #   vaccine_safety_claim OR treatment_efficacy_claim
        # ("system", "..."),
        ("user", "{post_text}"),
    ])
    # classifier_chain = ...  (compose prompt with llm)

    # (b) Router agent with tools. Must choose EXACTLY ONE tool and pass BOTH args.
    router_tools = [route_to_vaccine_fact_check, route_to_treatment_fact_check]
    for t in router_tools:
        t.handle_tool_error = True  # keep; agent should self-correct on bad calls

    router_prompt = ChatPromptTemplate.from_messages([
        # ("system", "... routing policy + pass both args exactly ..."),
        ("human",
         "Routing inputs:\n"
         "claim_type: {claim_type}\n"
         "is_priority: {is_priority}\n"
         "post: {post}"),
        MessagesPlaceholder("agent_scratchpad"),
    ])
    # router_agent = ...
    # router_executor = ...

    if USE_BLUESKY_DATA:
        print(f"\n⚠️ Live mode. Engagement threshold: {ENGAGEMENT_THRESHOLD}")

    # -------- Process each post --------
    for i, post in enumerate(posts):
        print("\n" + "-"*60)
        print(f"Processing Post {i+1}/{len(posts)} (ID: {post['id']})")
        print(f'Text: "{post["text"]}"')
        print(f"Engagement: {post['engagement']}")

        # -------- STUDENT TODO 2: PRIORITY + CLASSIFICATION --------
        # is_priority = ...
        # result = classifier_chain.invoke({...})
        # raw = ...
        # claim_type = ...  # normalize to exactly one of the two labels

        # print(f"  🚦 Priority: {'HIGH' if is_priority else 'Normal'}  |  🧠 Classified: {claim_type}")

        # -------- STUDENT TODO 3: ROUTE VIA AGENT --------
        # router_executor.invoke({
        #     "input": "Route now.",
        #     "post": post,
        #     "is_priority": is_priority,
        #     "claim_type": claim_type,
        # })

        time.sleep(0.2)

    fact_check_queues.display_queues()

if __name__ == "__main__":
    main()

