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

Output ONLY the following format in English (do not use any Chinese characters):

SEARCH STRATEGY:
[PubMed search query with MeSH terms and field tags]

EXPLANATION:
[2-3 sentence explanation]"""
    
    completion = client.chat.completions.create(
        model="qwen-qwq-32b",
        messages=[
            {
                "role": "system",
                "content": "You are a medical librarian who creates precise PubMed search strategies. Output only the search strategy and brief explanation in English, no thinking process or Chinese characters."
            },
            {
                "role": "user",
                "content": prompt
            }
        ],
        temperature=0.2,
        max_tokens=10000  # Increased from 500
    )
    
    response = completion.choices[0].message.content
    
    # Remove any thinking process enclosed in <think> tags and extra spaces
    response = re.sub(r'<think>.*?</think>\s*', '', response, flags=re.DOTALL)
    # Remove any leading/trailing whitespace
    response = response.strip()
    
    return response

def generate_summary(papers, research_question=None, max_batch_size=10):
    """
    Generate a holistic summary of the papers that addresses the research question.
    
    Parameters:
    -----------
    papers : list
        List of dictionaries containing paper information
    research_question : str, optional
        The research question being investigated
    max_batch_size : int, optional
        Maximum number of papers to process in one batch (default: 10)
        
    Returns:
    --------
    str
        A structured summary of the literature
    """
    if not papers:
        return "No papers provided for summarization."
    
    def create_summary_prompt(batch_papers):
        # Prepare the papers data for the prompt
        papers_text = "\n\n".join([
            f"Title: {paper.get('title', 'No title')}\n"
            f"Abstract: {paper.get('abstract', 'No abstract')}"
            for paper in batch_papers
        ])
        
        return f"""Research Question: {research_question if research_question else 'Not provided'}

Papers to analyze:
{papers_text}

Based on the above papers, provide a comprehensive analysis with the following structure. Use English only, do not output any Chinese characters:

Overall Summary:
- Focus on answering the research question directly. If few or no studies answer the research question, state this clearly and provide a very brief overview of what the included studies have investigated.
- State the level of evidence as one of: "Weak", "Moderate", or "Strong". Elaborate briefly on the type of evidence (e.g., meta-analysis, RCTs, case reports).

Gaps in Knowledge:
- Focus specifically on gaps related to the research question.
- If few or no studies directly answer the research question, explicitly state this.
- Be concise and specific about what aspects need further investigation.

Summary:
- Provide a synthesis of all knowledge relevant to the research question.
- If there is insufficient data to support any conclusions, state this clearly.
- Focus on concrete findings and their implications.

Please ensure each section directly relates to the research question. Keep the response focused and concise.

If you need to think through your analysis, enclose your thinking process in <think></think> tags.
"""

    try:
        client = get_groq_client()
        
        # Try with initial batch size
        current_batch_size = max_batch_size
        success = False
        
        while not success and current_batch_size > 0:
            try:
                # Take the first n papers
                batch_papers = papers[:current_batch_size]
                context = create_summary_prompt(batch_papers)
                
                chat_completion = client.chat.completions.create(
                    messages=[
                        {
                            "role": "system",
                            "content": "You are a scientific literature analyst. Provide clear, concise, and structured summaries that focus on addressing the research question. Use English only, no Chinese characters. Always maintain the exact structure provided in the prompt."
                        },
                        {
                            "role": "user",
                            "content": context
                        }
                    ],
                    model="qwen-qwq-32b",
                    max_tokens=10000,  
                    temperature=0.2,
                )
                
                summary = chat_completion.choices[0].message.content
                
                # Remove any thinking process enclosed in <think> tags and extra spaces
                summary = re.sub(r'<think>.*?</think>\s*', '', summary, flags=re.DOTALL)
                # Remove any leading/trailing whitespace
                summary = summary.strip()
                
                # If there are more papers that weren't included in this batch
                if len(papers) > current_batch_size:
                    summary += f"\n\nNote: This summary is based on {current_batch_size} most recent papers out of {len(papers)} total papers."
                
                success = True
                return summary
                
            except Exception as e:
                if "context window" in str(e).lower() or "api" in str(e).lower():
                    # Reduce batch size and try again
                    current_batch_size -= 1
                else:
                    # If it's a different type of error, raise it
                    raise e
        
        if not success:
            return "Unable to generate summary. The papers may be too long for the model's context window."
        
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

Analyze these {len(papers_to_analyze)} most recent papers to answer the question. Use English only, do not output any Chinese characters:
{papers_text}

Provide a concise response in this format:
ANSWER: [Brief, evidence-based answer focusing only on the question asked]
EVIDENCE: [Key findings from most relevant papers]
LIMITATIONS: [Main limitation in 1-2 sentences]

If you need to think through your analysis, enclose your thinking process in <think></think> tags.
"""

    try:
        client = get_groq_client()
        
        chat_completion = client.chat.completions.create(
            messages=[
                {
                    "role": "system",
                    "content": "You are a focused scientific analyst. Provide brief, evidence-based answers in English only, no Chinese characters."
                },
                {
                    "role": "user",
                    "content": context
                }
            ],
            model="qwen-qwq-32b",
            max_tokens=20000, 
            temperature=0.2,
        )
        
        response = chat_completion.choices[0].message.content
        
        # Remove any thinking process enclosed in <think> tags and extra spaces
        response = re.sub(r'<think>.*?</think>\s*', '', response, flags=re.DOTALL)
        # Remove any leading/trailing whitespace
        response = response.strip()
        
        # Add note about limited analysis
        if len(papers) > max_papers:
            response += f"\n\nNote: This analysis is based on the {max_papers} most recent papers out of {len(papers)} total papers."
        
        return response
        
    except Exception as e:
        raise Exception(f"Error analyzing literature: {str(e)}") 