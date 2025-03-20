import os
import re
import streamlit as st
from groq import Groq
from dotenv import load_dotenv

# Load environment variables
load_dotenv()

def get_groq_client():
    """Initialize and return a Groq client."""
    api_key = st.session_state.get('groq_api_key') or os.getenv('GROQ_API_KEY')
    if not api_key:
        raise ValueError("GROQ_API_KEY not found. Please enter your API key in the field above.")
    return Groq(api_key=api_key)

def generate_search_strategy(research_question: str) -> str:
    """Generate a PubMed search strategy using the QwQ-32B model."""
    client = get_groq_client()
    
    prompt = f"""Create a PubMed search strategy for: "{research_question}"

Output ONLY the following format, nothing else:

SEARCH STRATEGY:
[PubMed search query with MeSH terms and field tags]

EXPLANATION:
[2-3 sentence explanation]"""
    
    completion = client.chat.completions.create(
        model="qwen-qwq-32b",
        messages=[
            {
                "role": "system",
                "content": "You are a medical librarian who creates precise PubMed search strategies. Output only the search strategy and brief explanation, no thinking process."
            },
            {
                "role": "user",
                "content": prompt
            }
        ],
        temperature=0.2,
        max_tokens=1500  # Increased from 500
    )
    
    response = completion.choices[0].message.content
    
    # Remove any thinking process enclosed in <think> tags
    response = re.sub(r'<think>.*?</think>', '', response, flags=re.DOTALL)
    
    return response

def generate_summary(papers, research_question=None):
    """
    Generate a holistic summary of the papers that addresses the research question.
    
    Parameters:
    -----------
    papers : list
        List of dictionaries containing paper information
    research_question : str, optional
        The research question being investigated
        
    Returns:
    --------
    str
        A structured summary of the literature
    """
    if not papers:
        return "No papers provided for summarization."
    
    # Prepare the papers data for the prompt
    papers_text = "\n\n".join([
        f"Title: {paper.get('title', 'No title')}\n"
        f"Abstract: {paper.get('abstract', 'No abstract')}"
        for paper in papers
    ])
    
    # Construct the prompt
    context = f"""Research Question: {research_question if research_question else 'Not provided'}

Papers to analyze:
{papers_text}

Based on the above papers, provide a comprehensive analysis with the following structure:

1. Literature Overview:
- Synthesize the main findings and themes from the literature
- Evaluate how well the papers address the research question
- Identify any consensus or contradictions in the findings

2. Knowledge Gaps:
- Identify areas that are not well addressed in the current literature
- Point out limitations in the existing research
- Highlight questions that remain unanswered

3. Conclusions:
- Summarize the key takeaways in relation to the research question
- Assess the strength of the evidence
- Suggest directions for future research

Please ensure each section directly relates to the research question. If the papers do not adequately address the research question, explicitly state this and explain why.
"""

    try:
        client = get_groq_client()
        
        chat_completion = client.chat.completions.create(
            messages=[
                {
                    "role": "system",
                    "content": "You are a scientific literature analyst. Provide clear, concise, and structured summaries that focus on addressing the research question."
                },
                {
                    "role": "user",
                    "content": context
                }
            ],
            model="mixtral-8x7b-32768",
            max_tokens=3000,
            temperature=0.3,
        )
        
        return chat_completion.choices[0].message.content
        
    except Exception as e:
        raise Exception(f"Error generating summary: {str(e)}")

def ask_literature(papers, question: str, max_papers: int = 10) -> str:
    """
    Ask a specific question about the provided papers.
    
    Parameters:
    -----------
    papers : list
        List of dictionaries containing paper information
    question : str
        The specific question to answer based on the papers
    max_papers : int, optional
        Maximum number of papers to analyze (default: 10)
        
    Returns:
    --------
    str
        A focused answer to the question based on the literature
    """
    if not papers:
        return "No papers provided for analysis."
    
    # Sort papers by year (newest first) and limit the number
    sorted_papers = sorted(papers, key=lambda x: x.get('year', '0'), reverse=True)
    papers_to_analyze = sorted_papers[:max_papers]
    
    # Function to truncate abstract while keeping important parts
    def truncate_abstract(abstract: str, max_length: int = 300) -> str:
        if len(abstract) <= max_length:
            return abstract
        # Try to cut at the last complete sentence within limit
        truncated = abstract[:max_length]
        last_period = truncated.rfind('.')
        if last_period > 0:
            return abstract[:last_period + 1]
        return truncated + "..."
    
    # Prepare the papers data with truncated abstracts
    papers_text = "\n\n".join([
        f"Title: {paper.get('title', 'No title')}\n"
        f"Year: {paper.get('year', 'N/A')}\n"
        f"Abstract: {truncate_abstract(paper.get('abstract', 'No abstract'))}"
        for paper in papers_to_analyze
    ])
    
    # More focused prompt
    context = f"""Question: {question}

Analyze these {len(papers_to_analyze)} most recent papers to answer the question:
{papers_text}

Provide a concise response in this format:
ANSWER: [Brief, evidence-based answer focusing only on the question asked]
EVIDENCE: [Key findings from 2-3 most relevant papers]
LIMITATIONS: [Main limitation in 1-2 sentences]"""

    try:
        client = get_groq_client()
        
        chat_completion = client.chat.completions.create(
            messages=[
                {
                    "role": "system",
                    "content": "You are a focused scientific analyst. Provide brief, evidence-based answers."
                },
                {
                    "role": "user",
                    "content": context
                }
            ],
            model="mixtral-8x7b-32768",
            max_tokens=1000,  # Reduced from 2000
            temperature=0.3,
        )
        
        response = chat_completion.choices[0].message.content
        
        # Add note about limited analysis
        if len(papers) > max_papers:
            response += f"\n\nNote: This analysis is based on the {max_papers} most recent papers out of {len(papers)} total papers."
        
        return response
        
    except Exception as e:
        raise Exception(f"Error analyzing literature: {str(e)}") 