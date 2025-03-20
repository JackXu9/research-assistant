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
from pubmed_utils import search_pubmed, fetch_details
from groq_utils import generate_search_strategy, generate_summary, ask_literature
from search_utils import SearchManager
import ssl_bypass

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
    
    # Clear specific text input keys
    st.session_state.research_question = ""
    st.session_state.search_query = ""
    st.session_state.literature_question = ""
    
    # Clear multiselect values
    st.session_state["article_types"] = []  # Clear article types selection
    st.session_state["languages"] = []      # Clear languages selection

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
            # Add API key input at the top of the app
            groq_api_key = st.text_input("Enter your Groq API Key:", type="password")
            if groq_api_key:
                st.session_state.groq_api_key = groq_api_key
                os.environ["GROQ_API_KEY"] = groq_api_key
            
            # Research Question Section
            st.header("1. Define Research Question")
            research_question = st.text_area(
                "Enter your research question:", 
                key="research_question",
                placeholder="Example: Is kidney size associated with kidney function?"
            )
            
            if research_question:
                if st.button("Generate Search Strategy"):
                    with st.spinner("Generating search strategy..."):
                        try:
                            strategy = generate_search_strategy(research_question)
                            st.session_state.search_manager.current_query = strategy
                            st.text_area("Generated Search Strategy:", strategy, height=200)
                        except Exception as e:
                            st.error(f"Error generating search strategy: {str(e)}")

            # PubMed Search Section
            st.header("2. PubMed Search")
            
            # Email input
            email = st.text_input("Enter your email (required to search PubMed):", key="email")
            
            # Search options with fixed date range
            col1, col2 = st.columns(2)
            with col1:
                date_from = st.date_input("Date From:", datetime(2015, 1, 1), min_value=datetime(1900, 1, 1), max_value=datetime.now())
            with col2:
                date_to = st.date_input("Date To:", datetime.now(), min_value=datetime(1900, 1, 1), max_value=datetime.now())
            
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
            
            # Search query input
            search_query = st.text_area("Enter your PubMed search query:", key="search_query")
            
            if st.button("Execute Search"):
                if not email:
                    st.error("Please enter your PubMed email")
                elif not search_query:
                    st.error("Please enter a search query")
                else:
                    with st.spinner("Searching PubMed..."):
                        try:
                            # Build the base query first
                            final_query = search_query.strip()
                            
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
                            
                            # Format date with proper syntax
                            if date_from and date_to:
                                date_filter = f'{date_from.strftime("%Y/%m/%d")}:{date_to.strftime("%Y/%m/%d")}[Date - Publication]'
                                filter_components.append(date_filter)
                            
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
                            
                            # Add detailed debug output
                            st.write("Debug Information:")
                            st.write("1. Original Base Query:", original_query)
                            st.write("2. Applied Filters:")
                            if selected_languages:
                                st.write("   - Language filter:", ' OR '.join(language_filters))
                            if date_from and date_to:
                                st.write("   - Date filter:", date_filter)
                            if selected_article_types:
                                st.write("   - Article type filter:", ' OR '.join(type_filters))
                            st.write("3. Final Query:", final_query)
                            
                            # Add a link to test the base query without filters
                            encoded_base_query = original_query.replace(' ', '+')
                            base_pubmed_url = f"https://pubmed.ncbi.nlm.nih.gov/?term={encoded_base_query}"
                            st.markdown("4. Debug Links:")
                            st.markdown(f"   - [Test base query without filters]({base_pubmed_url})")
                            
                            # Execute search with correct parameters
                            results, count = search_pubmed(final_query, email, retmax=1000, verify_ssl=False)
                            
                            if results:
                                # Fetch details for each result
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
                                    
                                    st.success(f"Found {len(papers['PubmedArticle'])} papers")
                                else:
                                    st.warning("No valid paper data returned from PubMed")
                            else:
                                st.warning(f"No results found for query: {final_query}")
                                encoded_query = final_query.replace(' ', '+')
                                pubmed_url = f"https://pubmed.ncbi.nlm.nih.gov/?term={encoded_query}"
                                st.markdown(f"[Test full query on PubMed]({pubmed_url})")
                                st.info("""Troubleshooting suggestions:
        1. Try the base query link above without filters first
        2. If the base query works, add filters one at a time
        3. Check if the date range is too restrictive
        4. Verify that the MeSH terms exist in PubMed's vocabulary""")
                        except Exception as e:
                            st.error(f"Error during search: {str(e)}\nQuery: {final_query if 'final_query' in locals() else 'not constructed'}")

            # Results and Filtering Section
            if st.session_state.unfiltered_df is not None:
                st.header("3. Filter Results")
                
                # Filter by year
                years = sorted(st.session_state.unfiltered_df['year'].unique())
                selected_years = st.multiselect("Filter by Year:", years)
                
                # Filter by journal
                journals = sorted(st.session_state.unfiltered_df['journal'].unique())
                selected_journals = st.multiselect("Filter by Journal:", journals)
                
                # Apply filters
                if selected_years or selected_journals:
                    filtered_df = st.session_state.unfiltered_df.copy()
                    
                    if selected_years:
                        filtered_df = filtered_df[filtered_df['year'].isin(selected_years)]
                    
                    if selected_journals:
                        filtered_df = filtered_df[filtered_df['journal'].isin(selected_journals)]
                    
                    st.session_state.filtered_df = filtered_df
                
                # Summary Generation (moved to left column)
                st.header("4. Generate Summary")
                
                # Summary buttons in two columns for better spacing
                sum_col1, sum_col2 = st.columns(2)
                with sum_col1:
                    if st.button("Generate Full Summary", use_container_width=True):
                        if st.session_state.filtered_df is not None and not st.session_state.filtered_df.empty:
                            with st.spinner("Generating summary..."):
                                try:
                                    papers_to_summarize = st.session_state.filtered_df.to_dict('records')
                                    summary = generate_summary(papers_to_summarize, research_question)
                                    st.session_state.summary = summary
                                except Exception as e:
                                    st.error(f"Error generating summary: {str(e)}")
                        else:
                            st.warning("No papers selected for summarization")
                
                with sum_col2:
                    if st.button("Summarize Latest 10", use_container_width=True):
                        if st.session_state.filtered_df is not None and not st.session_state.filtered_df.empty:
                            with st.spinner("Generating summary of latest 10 papers..."):
                                try:
                                    latest_df = st.session_state.filtered_df.copy()
                                    latest_df['year'] = pd.to_numeric(latest_df['year'], errors='coerce')
                                    latest_df = latest_df.sort_values('year', ascending=False).head(10)
                                    
                                    if not latest_df.empty:
                                        papers_to_summarize = latest_df.to_dict('records')
                                        summary = generate_summary(papers_to_summarize, research_question)
                                        st.session_state.summary = summary
                                    else:
                                        st.warning("No valid papers found for summarization")
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
                    st.text_area("Generated Summary:", st.session_state.summary, height=400)
                
                # Ask the Literature section
                st.header("5. Ask the Literature")
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
                            st.markdown(f"<h3 class='study-title'>{row['title']}</h3>", unsafe_allow_html=True)
                            with st.container():
                                col1, col2 = st.columns([3, 1])
                                with col1:
                                    st.write("**Journal:**", row['journal'])
                                with col2:
                                    st.write("**PMID:**", row['pmid'])
                                
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
                    
                    # Create expandable section for each year
                    with st.expander(f"📅 {int(year)} ({len(year_papers)} papers)", expanded=True):
                        # Display each paper within the year group using a cleaner format
                        for _, row in year_papers.iterrows():
                            st.markdown(f"<h3 class='study-title'>{row['title']}</h3>", unsafe_allow_html=True)
                            with st.container():
                                col1, col2 = st.columns([3, 1])
                                with col1:
                                    st.write("**Journal:**", row['journal'])
                                with col2:
                                    st.write("**PMID:**", row['pmid'])
                                
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

if __name__ == "__main__":
    main() 