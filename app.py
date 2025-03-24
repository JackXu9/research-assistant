"""
Your Research Assistant

This application allows users to search PubMed for scientific papers
and generate summaries of the search results using the QwQ-32B reasoning model.
"""

import streamlit as st
import pandas as pd
from datetime import datetime
import os
from dotenv import load_dotenv
from pubmed_utils import search_pubmed, fetch_details, convert_ovid_to_pubmed
from groq_utils import generate_search_strategy, generate_summary, ask_literature, get_groq_client, generate_and_validate_search_strategy
from search_utils import SearchManager
import ssl_bypass
import re
import traceback

# Load environment variables
load_dotenv()

# Initialize session state
if 'search_manager' not in st.session_state:
    st.session_state.search_manager = SearchManager()
if 'unfiltered_df' not in st.session_state:
    st.session_state.unfiltered_df = None
if 'filtered_df' not in st.session_state:
    st.session_state.filtered_df = None
if 'summary' not in st.session_state:
    st.session_state.summary = None
if 'excluded_pmids' not in st.session_state:
    st.session_state.excluded_pmids = set()
if 'search_query' not in st.session_state:
    st.session_state.search_query = ""
# Initialize AI filter related session state
if 'ai_filtered_pmids' not in st.session_state:
    st.session_state.ai_filtered_pmids = set()
if 'ai_filter_explanations' not in st.session_state:
    st.session_state.ai_filter_explanations = {}

def reset_app():
    """Reset all session state variables except API key and email"""
    # Store current API key and email
    api_key = st.session_state.get('groq_api_key', '')
    email = st.session_state.get('email', '')
    
    # Clear all session state
    for key in list(st.session_state.keys()):
        del st.session_state[key]
    
    # Restore API key and email
    if api_key:
        st.session_state.groq_api_key = api_key
    if email:
        st.session_state.email = email
    
    # Reinitialize necessary session state variables
    st.session_state.search_manager = SearchManager()
    st.session_state.unfiltered_df = None
    st.session_state.filtered_df = None
    st.session_state.summary = None
    st.session_state.excluded_pmids = set()
    st.session_state.ai_filtered_pmids = set()
    st.session_state.ai_filter_explanations = {}
    
    # Clear specific text input keys
    st.session_state.research_question = ""
    st.session_state.search_query = ""
    st.session_state.literature_question = ""
    
    # Clear multiselect values
    st.session_state["article_types"] = []  # Clear article types selection
    st.session_state["languages"] = []      # Clear languages selection

def update_search_query(query):
    """Update the search query and trigger rerun"""
    st.session_state.search_query = query
    st.rerun()

# Move this BEFORE the main() function, right after your imports
st.set_page_config(
    page_title="Your Research Assistant",
    page_icon="📚",
    menu_items={
        'Get Help': None,
        'Report a bug': None,
        'About': """
# PubMed Summary App
Search PubMed and generate summaries using the QwQ-32B reasoning model.

---
**To reset the app**: 
1. Click the "⋮" menu in the top-right
2. Select "Rerun"
This will clear all inputs except API key and email.
        """
    }
)

def main():
    # Add sidebar information
    with st.sidebar:
        st.markdown("## How to get Groq API key")
        st.markdown("""
        Create a free Groq Cloud API-key at [console.groq.com/keys](https://console.groq.com/keys). 
        This allows for summarization and questioning of up to ten studies. 
        Processing of additional studies within the same batch can be achieved 
        with a pay per use subscription to groq cloud.
        """)
        
        # Add separator before reset button
        st.markdown("---")
        
        # Add Reset button
        if st.button("🔄 Reset App", key="reset_button"):
            reset_app()
            st.rerun()

    # Initialize session state if needed
    if 'initialized' not in st.session_state:
        reset_app()
        st.session_state.initialized = True

    # Add custom CSS to control layout width and positioning
    st.markdown("""
        <style>
        .main > div {
            max-width: 1500px;
            padding-left: 2rem;
            padding-right: 2rem;
            padding-top: 0.5rem;
        }
        .stApp {
            max-width: 100%;
        }
        .block-container {
            max-width: 1500px;
            padding-left: 2rem;
            padding-right: 2rem;
            padding-top: 1rem;
        }
        .element-container, .stTextArea > div > div > textarea {
            max-width: 100% !important;
        }
        .stButton > button {
            width: 100%;
        }
        .footer {
            position: fixed;
            bottom: 0;
            width: 100%;
            padding: 10px;
            text-align: center;
            font-size: 0.8em;
            border-top: 1px solid var(--st-color-border-light);
            background-color: var(--st-color-background);
            color: var(--st-color-text);
        }
        /* Reduce study title font size */
        .study-title {
            font-size: 1.17em !important;  /* Reduced from default h3 size */
        }
        </style>
    """, unsafe_allow_html=True)

    # Main app content
    st.title("Your Research Assistant")
    st.write("Search PubMed and generate summaries using the QwQ-32B reasoning model")

    # Create main layout columns with fixed width
    left_col, right_col = st.columns([0.6, 0.4])  # Adjusted back to original ratio

    with left_col:
        # Add container to control width
        with st.container():
            # Add API key and email inputs at the top of the app
            col1, col2 = st.columns(2)
            with col1:
                groq_api_key = st.text_input("Enter your Groq API Key:", type="password")
                if groq_api_key:
                    st.session_state.groq_api_key = groq_api_key
                    os.environ["GROQ_API_KEY"] = groq_api_key
            
            with col2:
                email = st.text_input("Enter your email (required for PubMed):", key="email")
            
            # Research Question Section
            st.header("1. Define Research Question")
            research_question = st.text_area(
                "Enter your research question:", 
                key="research_question",
                placeholder="Example: Is kidney size associated with kidney function?"
            )
            
            if research_question:
                if st.button("Search PubMed", type="primary"):
                    if not email:
                        st.error("Please enter your email address (required for PubMed)")
                    elif not st.session_state.get('groq_api_key'):
                        st.error("Please enter your Groq API key")
                    else:
                        # Start the auto-search process
                        if 'auto_search_state' not in st.session_state:
                            st.session_state.auto_search_state = {
                                "current_iteration": 0,
                                "max_iterations": 5,
                                "original_query": None,
                                "current_query": None,
                                "results": None,
                                "count": 0,
                                "all_attempts": [],
                                "success": False,
                                "error": None
                            }
                        
                        auto_search_progress = st.progress(0)
                        
                        with st.spinner("Generating search strategy and searching PubMed..."):
                            try:
                                # Reset the state for a new search
                                st.session_state.auto_search_state = {
                                    "current_iteration": 1,
                                    "max_iterations": 5,
                                    "original_query": None,
                                    "current_query": None,
                                    "results": None,
                                    "count": 0,
                                    "all_attempts": [],
                                    "success": False,
                                    "error": None
                                }
                                
                                # Generate initial search strategy
                                validated_strategy, validation_notes, was_modified = debug_generate_and_validate_search_strategy(research_question, email)
                                
                                # Extract the search query using regex
                                try:
                                    import re  # Explicitly import re here
                                    pattern = r'SEARCH STRATEGY:\s*(.*?)(?:\s*EXPLANATION:|$)'
                                    search_query_match = re.search(pattern, validated_strategy, re.DOTALL)
                                except Exception as e:
                                    st.error(f"Error in regex extraction: {str(e)}")
                                    st.write(f"Debug: Traceback: {traceback.format_exc()}")
                                    search_query_match = None
                                
                                if search_query_match:
                                    query = search_query_match.group(1).strip()
                                    st.session_state.auto_search_state["original_query"] = query
                                    st.session_state.auto_search_state["current_query"] = query
                                    
                                    # Save this attempt
                                    st.session_state.auto_search_state["all_attempts"].append({
                                        "iteration": 1,
                                        "query": query,
                                        "notes": validation_notes if was_modified else ["Initial query generation"]
                                    })
                                    
                                    st.info(f"Generated search query:")
                                    st.code(query)
                                    
                                    auto_search_progress.progress(0.1)
                                    
                                    # Now execute the search
                                    execute_auto_search(query, email, auto_search_progress)
                                else:
                                    st.error("Could not extract search query from the generated strategy")
                            except Exception as e:
                                st.error(f"Error in search process: {str(e)}")
                        
                        # Show results if the search was successful
                        if st.session_state.auto_search_state["success"]:
                            st.success(f"Search successful! Found {st.session_state.auto_search_state['count']} results.")
                            st.session_state.search_query = st.session_state.auto_search_state["current_query"]
                            
                            # Add Widen/Narrow search buttons
                            refine_col1, refine_col2 = st.columns(2)
                            with refine_col1:
                                if st.button("Widen Search", key="widen_search_btn"):
                                    if not email:
                                        st.error("Please enter your email address (required for PubMed)")
                                    elif not st.session_state.get('groq_api_key'):
                                        st.error("Please enter your Groq API key")
                                    else:
                                        with st.spinner("Widening search criteria..."):
                                            try:
                                                # Initialize Groq client
                                                client = get_groq_client()
                                                
                                                # Current query and research question
                                                current_query = st.session_state.auto_search_state["current_query"]
                                                current_count = st.session_state.auto_search_state["count"]
                                                
                                                # Create a prompt for widening search
                                                widen_prompt = f"""You are an expert PubMed search specialist. 
A search with the following query returned only {current_count} results, which is too few:

RESEARCH QUESTION: {research_question}

CURRENT QUERY:
{current_query}

Please create a more general, broader version of this query that will return MORE results.
Make these changes:
1. Use more general MeSH terms (move up the MeSH hierarchy)
2. Add synonyms with OR operators
3. Remove restrictive filters or qualifiers
4. Include related concepts
5. Ensure proper PubMed syntax

Return ONLY the revised search query with no explanation or additional text."""
                                                
                                                # Get the response from Groq
                                                completion = client.chat.completions.create(
                                                    messages=[
                                                        {
                                                            "role": "system",
                                                            "content": "You are a PubMed search expert who creates valid search queries. Respond with just the refined query, no explanation or preamble."
                                                        },
                                                        {
                                                            "role": "user",
                                                            "content": widen_prompt
                                                        }
                                                    ],
                                                    model="qwen-qwq-32b",
                                                    temperature=0.2,
                                                )
                                                
                                                # Extract the refined query
                                                import re
                                                widened_query = completion.choices[0].message.content.strip()
                                                
                                                # Clean up any markdown code blocks
                                                if "```" in widened_query:
                                                    try:
                                                        code_block_match = re.search(r'```(?:\w+\n)?(.*?)```', widened_query, re.DOTALL)
                                                        if code_block_match:
                                                            widened_query = code_block_match.group(1).strip()
                                                    except Exception:
                                                        # If regex fails, do basic string cleanup
                                                        widened_query = widened_query.replace("```", "").strip()
                                                
                                                # Reset the search state
                                                st.session_state.auto_search_state = {
                                                    "current_iteration": 1,
                                                    "max_iterations": 5,
                                                    "original_query": widened_query,
                                                    "current_query": widened_query,
                                                    "results": None,
                                                    "count": 0,
                                                    "all_attempts": [
                                                        {
                                                            "iteration": 1,
                                                            "query": widened_query,
                                                            "notes": ["Widened search criteria"]
                                                        }
                                                    ],
                                                    "success": False,
                                                    "error": None
                                                }
                                                
                                                # Display the new query
                                                st.info("Generated a wider search query:")
                                                st.code(widened_query)
                                                
                                                # Execute the new search
                                                execute_auto_search(widened_query, email)
                                                
                                                # Force a rerun to update the interface
                                                st.rerun()
                                                
                                            except Exception as e:
                                                st.error(f"Error widening search: {str(e)}")
                            
                            with refine_col2:
                                if st.button("Narrow Search", key="narrow_search_btn"):
                                    if not email:
                                        st.error("Please enter your email address (required for PubMed)")
                                    elif not st.session_state.get('groq_api_key'):
                                        st.error("Please enter your Groq API key")
                                    else:
                                        with st.spinner("Narrowing search criteria..."):
                                            try:
                                                # Initialize Groq client
                                                client = get_groq_client()
                                                
                                                # Current query and research question
                                                current_query = st.session_state.auto_search_state["current_query"]
                                                current_count = st.session_state.auto_search_state["count"]
                                                
                                                # Create a prompt for narrowing search
                                                narrow_prompt = f"""You are an expert PubMed search specialist. 
A search with the following query returned {current_count} results, which is too many:

RESEARCH QUESTION: {research_question}

CURRENT QUERY:
{current_query}

Please create a more specific, focused version of this query that will return FEWER but more relevant results.
Make these changes:
1. Use more specific MeSH terms (move down the MeSH hierarchy)
2. Add more specific qualifiers
3. Add key filters (publication types, dates, languages)
4. Focus on the core concepts of the research question 
5. Ensure proper PubMed syntax

Return ONLY the revised search query with no explanation or additional text."""
                                                
                                                # Get the response from Groq
                                                completion = client.chat.completions.create(
                                                    messages=[
                                                        {
                                                            "role": "system",
                                                            "content": "You are a PubMed search expert who creates valid search queries. Respond with just the refined query, no explanation or preamble."
                                                        },
                                                        {
                                                            "role": "user",
                                                            "content": narrow_prompt
                                                        }
                                                    ],
                                                    model="qwen-qwq-32b",
                                                    temperature=0.2,
                                                )
                                                
                                                # Extract the refined query
                                                import re
                                                narrowed_query = completion.choices[0].message.content.strip()
                                                
                                                # Clean up any markdown code blocks
                                                if "```" in narrowed_query:
                                                    try:
                                                        code_block_match = re.search(r'```(?:\w+\n)?(.*?)```', narrowed_query, re.DOTALL)
                                                        if code_block_match:
                                                            narrowed_query = code_block_match.group(1).strip()
                                                    except Exception:
                                                        # If regex fails, do basic string cleanup
                                                        narrowed_query = narrowed_query.replace("```", "").strip()
                                                
                                                # Reset the search state
                                                st.session_state.auto_search_state = {
                                                    "current_iteration": 1,
                                                    "max_iterations": 5,
                                                    "original_query": narrowed_query,
                                                    "current_query": narrowed_query,
                                                    "results": None,
                                                    "count": 0,
                                                    "all_attempts": [
                                                        {
                                                            "iteration": 1,
                                                            "query": narrowed_query,
                                                            "notes": ["Narrowed search criteria"]
                                                        }
                                                    ],
                                                    "success": False,
                                                    "error": None
                                                }
                                                
                                                # Display the new query
                                                st.info("Generated a narrower search query:")
                                                st.code(narrowed_query)
                                                
                                                # Execute the new search
                                                execute_auto_search(narrowed_query, email)
                                                
                                                # Force a rerun to update the interface
                                                st.rerun()
                                                
                                            except Exception as e:
                                                st.error(f"Error narrowing search: {str(e)}")
                            
                            # Display a summary of the search
                            with st.expander("Search Details", expanded=True):
                                for i, attempt in enumerate(st.session_state.auto_search_state["all_attempts"]):
                                    status = "✅ Success" if i == len(st.session_state.auto_search_state["all_attempts"]) - 1 and st.session_state.auto_search_state["success"] else "❌ Failed"
                                    st.markdown(f"**Attempt {attempt['iteration']}:** {status}")
                                    st.code(attempt["query"])
                                    if "result_count" in attempt:
                                        st.info(f"Results: {attempt['result_count']}")
                                    if "error" in attempt:
                                        st.error(f"Error: {attempt['error']}")
                                    if i < len(st.session_state.auto_search_state["all_attempts"]) - 1:
                                        st.markdown("---")
                        
                        # If we reached max iterations without success, show the error
                        elif st.session_state.auto_search_state["current_iteration"] >= st.session_state.auto_search_state["max_iterations"]:
                            st.error(f"Maximum {st.session_state.auto_search_state['max_iterations']} search attempts reached without success")
                            
                            # Show the last error if any
                            if st.session_state.auto_search_state["error"]:
                                st.error(f"Last error: {st.session_state.auto_search_state['error']}")
                            
                            # Display all attempts
                            with st.expander("All Search Attempts", expanded=True):
                                for attempt in st.session_state.auto_search_state["all_attempts"]:
                                    st.markdown(f"**Attempt {attempt['iteration']}:**")
                                    st.code(attempt["query"])
                                    if "result_count" in attempt:
                                        st.info(f"Results: {attempt['result_count']}")
                                    if "error" in attempt:
                                        st.error(f"Error: {attempt['error']}")
                                    st.markdown("---")

            # PubMed Search Section - changed to Filtering Results
            st.header("2. Filtering Results")
            
            if st.session_state.unfiltered_df is not None:
                # Show count of articles found
                st.info(f"Found {len(st.session_state.unfiltered_df)} articles. Apply filters below if needed.")
                
                # Filter by year
                years = sorted(st.session_state.unfiltered_df['year'].unique())
                selected_years = st.multiselect("Filter by Year:", years)
                
                # Article types (multiple selection)
                article_types = [
                    "Clinical Trial", "Randomized Controlled Trial", "Review", 
                    "Systematic Review", "Meta-Analysis", "Case Reports",
                    "Observational Study", "Cohort Study", "Cross-Sectional Study"
                ]
                selected_article_types = st.multiselect(
                    "Article Types:", 
                    article_types,
                    key="article_types"
                )
                
                # Languages (multiple selection)
                languages = ["eng", "fre", "ger", "spa", "ita", "por", "rus", "chi", "jpn", "kor"]
                language_names = {
                    "eng": "English", "fre": "French", "ger": "German", "spa": "Spanish",
                    "ita": "Italian", "por": "Portuguese", "rus": "Russian", "chi": "Chinese",
                    "jpn": "Japanese", "kor": "Korean"
                }
                selected_languages = st.multiselect(
                    "Languages:",
                    options=languages,
                    format_func=lambda x: language_names[x],
                    key="languages"
                )
                
                if st.button("Apply Filters"):
                    try:
                        # Apply DataFrame-based filters first (year)
                        if selected_years:
                            filtered_df = st.session_state.unfiltered_df.copy()
                            
                            if selected_years:
                                filtered_df = filtered_df[filtered_df['year'].isin(selected_years)]
                            
                            st.session_state.filtered_df = filtered_df
                            st.success(f"Filtered to {len(filtered_df)} papers based on year criteria")
                        
                        # Check if we need to apply PubMed-based filters (article types, languages)
                        if selected_languages or selected_article_types:
                            # Build the base query first
                            final_query = st.session_state.auto_search_state["current_query"]
                            
                            # Store original query for debugging
                            original_query = final_query
                            
                            # Add filters one by one, with proper formatting
                            filter_components = []
                            
                            if selected_languages:
                                language_filters = [f'{lang}[Language]' for lang in selected_languages]
                                if len(language_filters) == 1:
                                    filter_components.append(language_filters[0])
                                else:
                                    filter_components.append(f"({' OR '.join(language_filters)})")
                            
                            if selected_article_types:
                                type_filters = [f'"{at}"[Publication Type]' for at in selected_article_types]
                                if len(type_filters) == 1:
                                    filter_components.append(type_filters[0])
                                else:
                                    filter_components.append(f"({' OR '.join(type_filters)})")
                            
                            # Add all filters to the query
                            if filter_components:
                                final_query = f"({final_query}) AND {' AND '.join(filter_components)}"
                            
                            # Remove any extra spaces and normalize
                            final_query = ' '.join(final_query.split())
                            
                            # Execute search with filters
                            with st.spinner('Searching PubMed with filters...'):
                                results, count, warnings, cleaned_query = search_pubmed(final_query, email, retmax=1000, verify_ssl=False)
                                
                                if warnings:
                                    with st.expander("⚠️ Search Query Adjustments Found", expanded=True):
                                        st.markdown("**Original Query:**")
                                        st.code(final_query)
                                        st.markdown("**Suggested Changes:**")
                                        for warning in warnings:
                                            st.info(warning)
                                        st.markdown("**Updated Query:**")
                                        st.code(cleaned_query)

                                if not results:
                                    st.warning(f"No results found with the applied filters. Try broadening your criteria.")
                                    return

                                # Process results if we have them
                                if results:
                                    # Fetch details for each result
                                    with st.spinner('Fetching paper details...'):
                                        papers = fetch_details(results, email, verify_ssl=False)
                                    
                                    if papers and 'PubmedArticle' in papers:
                                        # Convert to DataFrame
                                        df = pd.DataFrame([{
                                            'title': paper['MedlineCitation']['Article'].get('ArticleTitle', 'Not available'),
                                            'journal': paper['MedlineCitation']['Article']['Journal'].get('Title', 'Not available'),
                                            'year': paper['MedlineCitation']['Article']['Journal']['JournalIssue']['PubDate'].get('Year', 'Not available'),
                                            'pmid': paper['MedlineCitation']['PMID'],
                                            'abstract': " ".join(paper['MedlineCitation']['Article'].get('Abstract', {}).get('AbstractText', ['Not available']))
                                        } for paper in papers['PubmedArticle']])
                                        
                                        # Apply DataFrame filters again to ensure consistency
                                        if selected_years:
                                            df = df[df['year'].isin(selected_years)]
                                        
                                        st.session_state.filtered_df = df.copy()
                                        
                                        st.success(f"Found {len(df)} papers after applying all filters")
                                    else:
                                        st.warning("No valid paper data returned from PubMed after applying filters")
                    except Exception as e:
                        st.error(f"Error applying filters: {str(e)}")
                        st.write(f"Debug: {traceback.format_exc()}")
            else:
                st.info("Execute a search first to see filtering options.")

            # Results and Filtering Section
            if st.session_state.unfiltered_df is not None:
                st.header("3. Generate Summary")
                
                # Get non-excluded papers
                included_df = st.session_state.filtered_df[~st.session_state.filtered_df['pmid'].isin(st.session_state.excluded_pmids)]
                
                if st.button("Generate Summary", use_container_width=True):
                    if not included_df.empty:
                        with st.spinner("Generating summary..."):
                            try:
                                papers_to_summarize = included_df.to_dict('records')
                                summary = generate_summary(papers_to_summarize, research_question)
                                st.session_state.summary = summary
                            except Exception as e:
                                st.error(f"Error generating summary: {str(e)}")
                    else:
                        st.warning("No papers selected for summarization")
                
                # Download buttons in two columns
                down_col1, down_col2 = st.columns(2)
                with down_col1:
                    if st.session_state.summary:
                        st.download_button(
                            "Download Summary",
                            st.session_state.summary,
                            file_name="summary.txt",
                            mime="text/plain",
                            use_container_width=True
                        )
                
                with down_col2:
                    if st.session_state.filtered_df is not None:
                        csv = st.session_state.filtered_df.to_csv(index=False)
                        st.download_button(
                            "Download Search Results",
                            csv,
                            file_name="search_results.csv",
                            mime="text/csv",
                            use_container_width=True
                        )
                
                # Display summary text area
                if st.session_state.summary:
                    # Format specific sections to be bold using regex
                    import re
                    formatted_summary = st.session_state.summary
                    
                    # Format the Level of Evidence section
                    formatted_summary = re.sub(
                        r'(Level of [Ee]vidence:?|[Ee]vidence [Ll]evel:?|[Ee]vidence [Qq]uality:?)(.*?)(\n|$)', 
                        r'**\1**\2\3', 
                        formatted_summary
                    )
                    
                    # Format the Implications or Conclusion section
                    formatted_summary = re.sub(
                        r'([Ii]mplications:?|[Cc]onclusions?:?|[Cc]linical [Ii]mplications:?)(.*?)(\n|$)', 
                        r'**\1**\2\3', 
                        formatted_summary
                    )
                    
                    # Display the formatted summary
                    st.markdown(formatted_summary)
                
                # Ask the Literature section
                st.header("4. Ask the Literature")
                st.write("Ask a specific question about the selected papers.")
                
                # Question input
                literature_question = st.text_area("Enter your question:", key="literature_question", 
                                                 placeholder="e.g., What are the most common side effects reported in these studies?")
                
                # Create columns for the buttons
                ask_col1, ask_col2 = st.columns([0.7, 0.3])
                
                with ask_col1:
                    if st.button("Ask Question", use_container_width=True):
                        if not literature_question:
                            st.warning("Please enter a question first.")
                        elif st.session_state.filtered_df is not None and not st.session_state.filtered_df.empty:
                            with st.spinner("Analyzing literature..."):
                                try:
                                    papers_to_analyze = st.session_state.filtered_df.to_dict('records')
                                    answer = ask_literature(papers_to_analyze, literature_question)
                                    st.session_state.literature_answer = answer
                                except Exception as e:
                                    st.error(f"Error analyzing literature: {str(e)}")
                        else:
                            st.warning("No papers selected for analysis")
                
                with ask_col2:
                    if 'literature_answer' in st.session_state and st.session_state.literature_answer:
                        st.download_button(
                            "Download Answer",
                            st.session_state.literature_answer,
                            file_name="literature_analysis.txt",
                            mime="text/plain",
                            use_container_width=True
                        )
                
                # Display the answer if available
                if 'literature_answer' in st.session_state and st.session_state.literature_answer:
                    st.markdown("### Analysis Result")
                    st.markdown(st.session_state.literature_answer)

    # Right column for displaying papers
    with right_col:
        if st.session_state.filtered_df is not None:
            st.header("Search Results")
            df_display = st.session_state.filtered_df.copy()
            
            # Display articles count
            total_articles = len(df_display)
            if st.session_state.ai_filtered_pmids:
                filtered_count = len(st.session_state.ai_filtered_pmids)
                st.info(f"Showing {filtered_count} of {total_articles} articles")
            else:
                st.info(f"Showing all {total_articles} articles")
            
            # Add Reset All Filters button - keep for when filters are applied
            if st.button("🔄 Reset All Filters", key="reset_all_filters_btn", use_container_width=True):
                # Clear AI filter
                st.session_state.ai_filtered_pmids = set()
                st.session_state.ai_filter_explanations = {}
                
                # Reset manual filters (year selections)
                # Clear all multiselect filters in the left column
                for key in list(st.session_state.keys()):
                    if key.startswith('Filter by') or key.endswith('multiselect'):
                        st.session_state[key] = []
                    
                # Reset the filtered_df to the original unfiltered state
                if st.session_state.unfiltered_df is not None:
                    st.session_state.filtered_df = st.session_state.unfiltered_df.copy()
                
                st.success("All filters have been reset!")
                st.rerun()
                
            # Add AI Filter section
            st.subheader("AI Filter")
            ai_filter_criteria = st.text_area(
                "Enter filter criteria for the AI to analyze abstracts:",
                placeholder="Example: studies including populations more than 20 participants\nor: randomized controlled trials with follow-up > 6 months",
                help="The AI will analyze each abstract based on your criteria and filter the results accordingly."
            )
            
            # Initialize AI filtered papers in session state if not exists
            if 'ai_filtered_pmids' not in st.session_state:
                st.session_state.ai_filtered_pmids = set()
            
            # Initialize AI filter explanations
            if 'ai_filter_explanations' not in st.session_state:
                st.session_state.ai_filter_explanations = {}
            
            if ai_filter_criteria and st.button("Apply AI Filter"):
                # Check if API key is available
                if not st.session_state.get('groq_api_key'):
                    st.error("Please enter your Groq API key at the top of the page to use the AI filter.")
                else:
                    with st.spinner("AI analyzing abstracts..."):
                        try:
                            # Initialize Groq client
                            client = get_groq_client()
                            
                            # Prepare a more structured filtering prompt that returns JSON
                            filter_prompt = f"""Analyze the following abstract to determine if it meets this criteria:
{ai_filter_criteria}

IMPORTANT: Return your analysis as a JSON object with this structure:
{{
  "meets_criteria": true/false,
  "reason": "Brief explanation of why it does or doesn't meet the criteria",
  "extracted_info": {{
    "participant_count": number or null,
    "study_type": "string",
    "follow_up_period": "string or null"
  }}
}}

If you can't determine whether the abstract meets the criteria, set meets_criteria to false.
If you can't extract specific information, use null for that field.
Only include factual information directly from the abstract.

Abstract:
"""
                            filtered_pmids = set()
                            filter_explanations = {}
                            progress_bar = st.progress(0)
                            total_papers = len(df_display)
                            
                            for idx, row in df_display.iterrows():
                                # Create completion for each abstract
                                completion = client.chat.completions.create(
                                    messages=[
                                        {
                                            "role": "system",
                                            "content": "You are a scientific abstract analyzer that extracts structured information and determines if studies meet specific criteria. Always respond with valid JSON."
                                        },
                                        {
                                            "role": "user",
                                            "content": filter_prompt + row['abstract']
                                        }
                                    ],
                                    model="qwen-qwq-32b",
                                    temperature=0.1,
                                )
                                
                                # Process the response
                                response_text = completion.choices[0].message.content.strip()
                                try:
                                    # Extract the JSON part if wrapped in markdown code blocks
                                    if "```json" in response_text:
                                        response_text = response_text.split("```json")[1].split("```")[0].strip()
                                    elif "```" in response_text:
                                        response_text = response_text.split("```")[1].split("```")[0].strip()
                                    
                                    import json
                                    result = json.loads(response_text)
                                    
                                    # Store the PMID if it meets criteria
                                    if result.get("meets_criteria", False):
                                        pmid = str(row['pmid'])
                                        filtered_pmids.add(pmid)
                                        filter_explanations[pmid] = result.get("reason", "Meets criteria")
                                    else:
                                        # Store explanation for why it was excluded
                                        filter_explanations[str(row['pmid'])] = result.get("reason", "Does not meet criteria")
                                        
                                except Exception as json_e:
                                    print(f"Error parsing JSON response for PMID {row['pmid']}: {str(json_e)}")
                                    print(f"Response text: {response_text}")
                                    # Fallback to looking for simpler yes/no in the response
                                    if "true" in response_text.lower() or '"meets_criteria": true' in response_text.lower():
                                        filtered_pmids.add(str(row['pmid']))
                                        filter_explanations[str(row['pmid'])] = "Meets criteria (parsed from text)"
                                
                                # Update progress bar
                                progress_bar.progress((idx + 1) / total_papers)
                            
                            # Store results in session state
                            st.session_state.ai_filtered_pmids = filtered_pmids
                            st.session_state.ai_filter_explanations = filter_explanations
                            
                            # Show results
                            if filtered_pmids:
                                st.success(f"Found {len(filtered_pmids)} papers matching your criteria.")
                            else:
                                st.warning("No papers found matching your criteria.")
                            
                        except Exception as e:
                            st.error(f"Error during AI filtering: {str(e)}")
            
            # Add button to clear AI filter
            if st.session_state.ai_filtered_pmids:
                if st.button("Clear AI Filter", key="clear_ai_filter_btn"):
                    st.session_state.ai_filtered_pmids = set()
                    st.session_state.ai_filter_explanations = {}
                    st.rerun()
            
            # Group papers by year
            if not df_display.empty:
                # Convert year to numeric, handling any non-numeric values
                df_display['year'] = pd.to_numeric(df_display['year'], errors='coerce')
                # Sort years in descending order and remove NaN values
                years = sorted([y for y in df_display['year'].unique() if pd.notna(y)], reverse=True)
                
                # Add a section for papers with unknown years
                unknown_year_papers = df_display[df_display['year'].isna()]
                if not unknown_year_papers.empty:
                    with st.expander(f"📅 Unknown Year ({len(unknown_year_papers)} papers)", expanded=True):
                        for _, row in unknown_year_papers.iterrows():
                            # Skip if AI filter is active and paper doesn't match
                            if st.session_state.ai_filtered_pmids and str(row['pmid']) not in st.session_state.ai_filtered_pmids:
                                continue
                                
                            st.markdown(f"<h3 class='study-title'>{row['title']}</h3>", unsafe_allow_html=True)
                            with st.container():
                                col1, col2, col3 = st.columns([2.5, 0.8, 0.7])
                                with col1:
                                    st.write("**Journal:**", row['journal'])
                                with col2:
                                    st.write("**PMID:**", row['pmid'])
                                with col3:
                                    # Add toggle button for excluding/including the paper
                                    pmid = str(row['pmid'])
                                    if pmid in st.session_state.excluded_pmids:
                                        if st.button("➕ Include", key=f"include_{pmid}", type="primary"):
                                            st.session_state.excluded_pmids.remove(pmid)
                                            st.rerun()
                                    else:
                                        if st.button("❌ Exclude", key=f"exclude_{pmid}"):
                                            st.session_state.excluded_pmids.add(pmid)
                                            st.rerun()
                                
                                # Show filter reason if available
                                if pmid in st.session_state.ai_filter_explanations:
                                    st.info(f"**Filter match reason:** {st.session_state.ai_filter_explanations[pmid]}")
                                
                                # Add a toggle for the abstract
                                if st.button(f"Toggle Abstract 🔍", key=f"toggle_{row['pmid']}"):
                                    st.session_state[f'show_abstract_{row["pmid"]}'] = \
                                        not st.session_state.get(f'show_abstract_{row["pmid"]}', False)
                                
                                # Show abstract if toggled
                                if st.session_state.get(f'show_abstract_{row["pmid"]}', False):
                                    st.markdown(f"**Abstract:**\n{row['abstract']}")
                                
                            st.markdown("---")  # Add separator between papers
                
                # Display papers grouped by year
                for year in years:
                    year_papers = df_display[df_display['year'] == year]
                    
                    # Count papers that match AI filter if active
                    if st.session_state.ai_filtered_pmids:
                        matching_papers = year_papers[year_papers['pmid'].astype(str).isin(st.session_state.ai_filtered_pmids)]
                        paper_count = len(matching_papers)
                    else:
                        paper_count = len(year_papers)
                    
                    # Create expandable section for each year
                    with st.expander(f"📅 {int(year)} ({paper_count} papers)", expanded=True):
                        # Display each paper within the year group using a cleaner format
                        for _, row in year_papers.iterrows():
                            # Skip if AI filter is active and paper doesn't match
                            if st.session_state.ai_filtered_pmids and str(row['pmid']) not in st.session_state.ai_filtered_pmids:
                                continue
                                
                            st.markdown(f"<h3 class='study-title'>{row['title']}</h3>", unsafe_allow_html=True)
                            with st.container():
                                col1, col2, col3 = st.columns([2.5, 0.8, 0.7])
                                with col1:
                                    st.write("**Journal:**", row['journal'])
                                with col2:
                                    st.write("**PMID:**", row['pmid'])
                                with col3:
                                    # Add toggle button for excluding/including the paper
                                    pmid = str(row['pmid'])
                                    if pmid in st.session_state.excluded_pmids:
                                        if st.button("➕ Include", key=f"include_{pmid}", type="primary"):
                                            st.session_state.excluded_pmids.remove(pmid)
                                            st.rerun()
                                    else:
                                        if st.button("❌ Exclude", key=f"exclude_{pmid}"):
                                            st.session_state.excluded_pmids.add(pmid)
                                            st.rerun()
                                
                                # Show filter reason if available
                                if pmid in st.session_state.ai_filter_explanations:
                                    st.info(f"**Filter match reason:** {st.session_state.ai_filter_explanations[pmid]}")
                                
                                # Add a toggle for the abstract
                                if st.button(f"Toggle Abstract 🔍", key=f"toggle_{row['pmid']}"):
                                    st.session_state[f'show_abstract_{row["pmid"]}'] = \
                                        not st.session_state.get(f'show_abstract_{row["pmid"]}', False)
                                
                                # Show abstract if toggled
                                if st.session_state.get(f'show_abstract_{row["pmid"]}', False):
                                    st.markdown(f"**Abstract:**\n{row['abstract']}")
                                
                            st.markdown("---")  # Add separator between papers

    # Add footer at the bottom of the page
    st.markdown("""
        <div class='footer'>
            Created by Jack Xu, MD & PhD | Contact: jack.junchi.xu@regionh.dk
        </div>
    """, unsafe_allow_html=True)

def reset_and_rerun():
    """Helper function to reset and rerun the app"""
    reset_app()
    st.rerun()

def execute_auto_search(query, email, progress_bar=None):
    """
    Execute a PubMed search and handle refinement if needed.
    
    Args:
        query: The search query to execute
        email: User's email for PubMed
        progress_bar: Optional Streamlit progress bar to update
    
    Returns:
        None - updates session state
    """
    try:
        # Execute the search
        results, count, warnings, cleaned_query = search_pubmed(query, email, retmax=1000, verify_ssl=False)
        
        # Update the current attempt with results
        current_attempt = st.session_state.auto_search_state["all_attempts"][-1]
        current_attempt["result_count"] = count
        current_attempt["warnings"] = warnings
        current_attempt["cleaned_query"] = cleaned_query
        
        # Update progress
        if progress_bar:
            progress_bar.progress(0.3 + (0.7 * st.session_state.auto_search_state["current_iteration"] / st.session_state.auto_search_state["max_iterations"]))
        
        # Check if we got any results
        if count > 0:
            # Success - we have results
            st.session_state.auto_search_state["success"] = True
            st.session_state.auto_search_state["results"] = results
            st.session_state.auto_search_state["count"] = count
            st.session_state.auto_search_state["current_query"] = cleaned_query
            
            # Fetch paper details to prepare for display
            with st.spinner('Fetching paper details...'):
                papers = fetch_details(results, email, verify_ssl=False)
            
            if papers and 'PubmedArticle' in papers:
                # Convert to DataFrame
                df = pd.DataFrame([{
                    'title': paper['MedlineCitation']['Article'].get('ArticleTitle', 'Not available'),
                    'journal': paper['MedlineCitation']['Article']['Journal'].get('Title', 'Not available'),
                    'year': paper['MedlineCitation']['Article']['Journal']['JournalIssue']['PubDate'].get('Year', 'Not available'),
                    'pmid': paper['MedlineCitation']['PMID'],
                    'abstract': " ".join(paper['MedlineCitation']['Article'].get('Abstract', {}).get('AbstractText', ['Not available']))
                } for paper in papers['PubmedArticle']])
                
                st.session_state.unfiltered_df = df
                st.session_state.filtered_df = df.copy()
                
                # Update UI - this needs to be handled by the caller with a rerun
                st.success(f"Found {count} papers matching your search criteria.")
            
            return True
        else:
            # No results - need to refine
            if st.session_state.auto_search_state["current_iteration"] < st.session_state.auto_search_state["max_iterations"]:
                # Increment iteration counter
                st.session_state.auto_search_state["current_iteration"] += 1
                
                # Generate a more broad query
                refine_query(
                    st.session_state.auto_search_state["current_query"], 
                    0,  # Zero results
                    email, 
                    progress_bar
                )
            else:
                st.warning(f"No results found after {st.session_state.auto_search_state['max_iterations']} attempts")
            
            return False
    
    except Exception as e:
        # Store the error
        st.session_state.auto_search_state["error"] = str(e)
        
        # Update the current attempt with the error
        current_attempt = st.session_state.auto_search_state["all_attempts"][-1]
        current_attempt["error"] = str(e)
        
        # If we haven't hit max iterations, try to refine the query
        if st.session_state.auto_search_state["current_iteration"] < st.session_state.auto_search_state["max_iterations"]:
            # Increment iteration counter
            st.session_state.auto_search_state["current_iteration"] += 1
            
            # Generate a more broad query based on the error
            refine_query(
                st.session_state.auto_search_state["current_query"], 
                -1,  # Error code
                email,
                progress_bar,
                error_message=str(e)
            )
        else:
            st.error(f"Maximum attempts reached. Last error: {str(e)}")
        
        return False

def refine_query(current_query, result_count, email, progress_bar=None, error_message=None):
    """
    Refine a search query based on results or errors.
    
    Args:
        current_query: The current search query
        result_count: Number of results from previous search (-1 for error)
        email: User's email for PubMed
        progress_bar: Optional Streamlit progress bar to update
        error_message: Optional error message if there was an error
    
    Returns:
        None - updates session state and triggers new search
    """
    try:
        # Import re module at the beginning of the function
        import re
        
        # Initialize Groq client
        client = get_groq_client()
        
        # Create a refinement prompt based on the outcome
        if result_count == 0:
            refinement_prompt = f"""You are an expert PubMed search specialist. 
A PubMed search with the following query returned ZERO results:

```
{current_query}
```

Create a more broad, simplified version of this query that is more likely to return results.
Make these specific changes:
1. Remove any overly specific filters or qualifiers
2. Use more general terms
3. Simplify the boolean structure
4. Ensure proper PubMed syntax

Return ONLY the refined search query with no explanation or additional text."""
        else:  # Error occurred
            refinement_prompt = f"""You are an expert PubMed search specialist. 
A PubMed search with the following query returned an ERROR:

```
{current_query}
```

Error message: {error_message}

Rewrite this query to fix syntax errors and ensure it will work properly in PubMed. 
Make these specific changes:
1. Fix any syntax errors
2. Use proper PubMed field tags
3. Ensure balanced parentheses and quotes
4. Simplify the query structure if too complex

Return ONLY the corrected search query with no explanation or additional text."""
        
        # Get the response
        completion = client.chat.completions.create(
            messages=[
                {
                    "role": "system",
                    "content": "You are a PubMed search expert who creates valid search queries. Respond with just the refined query, no explanation or preamble."
                },
                {
                    "role": "user",
                    "content": refinement_prompt
                }
            ],
            model="qwen-qwq-32b",
            temperature=0.2,
        )
        
        # Extract the refined query
        refined_query = completion.choices[0].message.content.strip()
        
        # Clean up any markdown code blocks
        if "```" in refined_query:
            try:
                code_block_match = re.search(r'```(?:\w+\n)?(.*?)```', refined_query, re.DOTALL)
                if code_block_match:
                    refined_query = code_block_match.group(1).strip()
            except Exception as e:
                st.error(f"Error extracting code block: {str(e)}")
                # If regex fails, do basic string cleanup
                refined_query = refined_query.replace("```", "").strip()
        
        # Update session state
        st.session_state.auto_search_state["current_query"] = refined_query
        
        # Save this attempt
        st.session_state.auto_search_state["all_attempts"].append({
            "iteration": st.session_state.auto_search_state["current_iteration"],
            "query": refined_query,
            "notes": [f"Refined query after {'error' if result_count == -1 else 'zero results'}"]
        })
        
        # Update progress bar
        if progress_bar:
            progress_bar.progress(0.3 + (0.7 * (st.session_state.auto_search_state["current_iteration"] - 1) / st.session_state.auto_search_state["max_iterations"]))
        
        # Execute the refined search
        st.info(f"Trying refined query (Attempt {st.session_state.auto_search_state['current_iteration']} of {st.session_state.auto_search_state['max_iterations']})")
        st.code(refined_query)
        
        # Execute the search with the refined query
        execute_auto_search(refined_query, email, progress_bar)
        
    except Exception as e:
        st.error(f"Error refining query: {str(e)}")
        st.write(f"Debug traceback: {traceback.format_exc()}")
        st.session_state.auto_search_state["error"] = str(e)
        
        # If we still haven't hit max iterations, try once more with a simpler approach
        if st.session_state.auto_search_state["current_iteration"] < st.session_state.auto_search_state["max_iterations"]:
            # Increment iteration counter
            st.session_state.auto_search_state["current_iteration"] += 1
            
            try:
                # Generate an emergency simplified query
                emergency_query = emergency_simplify_query(current_query)
                
                # Update session state
                st.session_state.auto_search_state["current_query"] = emergency_query
                
                # Save this emergency attempt
                st.session_state.auto_search_state["all_attempts"].append({
                    "iteration": st.session_state.auto_search_state["current_iteration"],
                    "query": emergency_query,
                    "notes": ["Emergency simplified query after refinement error"]
                })
                
                # Execute the emergency search
                st.warning(f"Trying emergency simplified query (Attempt {st.session_state.auto_search_state['current_iteration']} of {st.session_state.auto_search_state['max_iterations']})")
                st.code(emergency_query)
                
                execute_auto_search(emergency_query, email, progress_bar)
            except Exception as e2:
                st.error(f"Emergency query also failed: {str(e2)}")
                st.write(f"Debug emergency traceback: {traceback.format_exc()}")

def emergency_simplify_query(query):
    """
    Create a drastically simplified version of a search query as a last resort.
    
    Args:
        query: The complex query that failed
    
    Returns:
        str: A basic simplified query
    """
    try:
        # Import re at function start
        import re
        
        # Extract main terms (words in quotes or words before field tags)
        terms = re.findall(r'"([^"]+)"', query)
        terms += re.findall(r'(\w+)\[', query)
        
        # Remove duplicates and filter out common words
        stopwords = ["and", "or", "not", "the", "a", "an", "in", "on", "of", "to", "for"]
        filtered_terms = [term for term in terms if len(term) > 2 and term.lower() not in stopwords]
        
        # Take up to 3 most likely meaningful terms
        main_terms = list(set(filtered_terms))[:3]
        
        # If we couldn't find any terms, extract words longer than 5 characters
        if not main_terms:
            words = re.findall(r'\b(\w{5,})\b', query)
            main_terms = list(set([w for w in words if w.lower() not in stopwords]))[:3]
        
        # If still empty, return a very basic query
        if not main_terms:
            return "heart disease[MeSH Terms]"
        
        # Create a simple OR query with the main terms as text words
        simple_query = " OR ".join([f'"{term}"[Text Word]' for term in main_terms])
        
        # Add a date range if we have one in the original query
        date_range = re.search(r'(\d{4}/\d{2}/\d{2}:\d{4}/\d{2}/\d{2})\[Date - Publication\]', query)
        if date_range:
            simple_query += f" AND {date_range.group(0)}"
        
        return simple_query
    except Exception as e:
        st.error(f"Error in emergency query simplification: {str(e)}")
        st.write(f"Debug traceback: {traceback.format_exc()}")
        # Return a very basic fallback query that should work
        return "heart[MeSH Terms]"

# Create a debug wrapper for the generate_and_validate_search_strategy function
def debug_generate_and_validate_search_strategy(research_question, email):
    """Debug wrapper for generate_and_validate_search_strategy with explicit error handling"""
    try:
        # Add explicit import of re inside the function
        import re
        return generate_and_validate_search_strategy(research_question, email)
    except Exception as e:
        st.error(f"Debug error in search strategy generation: {str(e)}")
        st.write(f"Traceback: {traceback.format_exc()}")
        # Return a fallback empty strategy
        return f"Error generating strategy: {str(e)}", ["Error occurred"], False

if __name__ == "__main__":
    main() 