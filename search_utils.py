"""
Search Utilities Module

This module provides utilities for building and managing PubMed search queries,
including a SearchManager class and functions for rendering search interfaces.
"""

import streamlit as st
import json
import os
from pathlib import Path

class SearchManager:
    """
    Class for managing and building PubMed search queries.
    
    This class maintains the current search query and provides methods
    for adding search terms with different operators.
    
    Attributes:
    -----------
    current_query : str
        The current search query string.
    """
    
    def __init__(self, initial_query=""):
        """Initialize the SearchManager with an optional initial query."""
        self.current_query = initial_query
    
    def add_term(self, term, field="", operator="AND"):
        """
        Add a term to the current query with the specified field and operator.
        
        Parameters:
        -----------
        term : str
            The search term to add.
        field : str, default=""
            The PubMed field to search in (e.g., "[Title/Abstract]").
        operator : str, default="AND"
            The Boolean operator to use (AND, OR, NOT).
        """
        if not term:
            return
            
        # Format the term with field if provided
        formatted_term = f"{term}{field}"
        
        # Add to the current query with the appropriate operator
        if not self.current_query:
            # First term doesn't need an operator
            self.current_query = formatted_term
        else:
            self.current_query = f"({self.current_query}) {operator} {formatted_term}"
    
    def reset_query(self):
        """Reset the current query to an empty string."""
        self.current_query = ""
    
    def save_query(self, name):
        """
        Save the current query with the given name.
        
        Parameters:
        -----------
        name : str
            The name to save the query under.
        """
        if not self.current_query or not name:
            return
            
        # Create data directory if it doesn't exist
        data_dir = Path("data")
        data_dir.mkdir(exist_ok=True)
        
        # Load existing saved searches
        saved_searches = {}
        saved_searches_path = data_dir / "saved_searches.json"
        if saved_searches_path.exists():
            try:
                with open(saved_searches_path, "r") as f:
                    saved_searches = json.load(f)
            except:
                saved_searches = {}
        
        # Add new search
        saved_searches[name] = self.current_query
        
        # Save updated searches
        with open(saved_searches_path, "w") as f:
            json.dump(saved_searches, f)

def render_advanced_search(search_manager):
    """
    Render the advanced search interface in Streamlit.
    
    Parameters:
    -----------
    search_manager : SearchManager
        The SearchManager instance to use for building the query.
    """
    st.write("Build your search by adding terms")
    
    cols = st.columns(4)
    with cols[0]:
        operator = st.selectbox("Operator", ["AND", "OR", "NOT"], index=0)
    with cols[1]:
        term = st.text_input("Term")
    with cols[2]:
        field = st.selectbox("Field", [
            "",
            "[Title/Abstract]",
            "[MeSH Terms]",
            "[Author]",
            "[Journal]",
            "[Publication Type]",
            "[Language]"
        ])
    with cols[3]:
        if st.button("Add Term"):
            search_manager.add_term(term, field, operator)
    
    # Date range filter
    st.write("### Date Range")
    date_cols = st.columns(2)
    with date_cols[0]:
        start_year = st.text_input("Start Year (YYYY)")
    with date_cols[1]:
        end_year = st.text_input("End Year (YYYY)")
    
    if st.button("Add Date Range") and (start_year or end_year):
        date_query = ""
        if start_year and end_year:
            date_query = f"{start_year}:{end_year}[dp]"
        elif start_year:
            date_query = f"{start_year}:3000[dp]"
        elif end_year:
            date_query = f"1800:{end_year}[dp]"
        
        if date_query:
            search_manager.add_term(date_query, "", "AND")
    
    # Language filter
    st.write("### Language")
    language = st.selectbox("Language", ["", "English", "French", "German", "Spanish", "Chinese", "Japanese"])
    if st.button("Add Language") and language:
        search_manager.add_term(f"{language}[Language]", "", "AND")
    
    # Article type filter
    st.write("### Article Type")
    article_type = st.selectbox("Article Type", [
        "",
        "Review",
        "Clinical Trial",
        "Meta-Analysis",
        "Randomized Controlled Trial",
        "Systematic Review",
        "Case Reports"
    ])
    if st.button("Add Article Type") and article_type:
        search_manager.add_term(f"{article_type}[Publication Type]", "", "AND")
    
    # Current query display
    st.write("### Current Query")
    st.code(search_manager.current_query)
    
    # Reset button
    if st.button("Reset Query"):
        search_manager.reset_query()
        st.rerun()
    
    # Save query
    st.write("### Save Query")
    query_name = st.text_input("Query Name")
    if st.button("Save") and query_name:
        search_manager.save_query(query_name)
        st.success(f"Query saved as '{query_name}'")
    
    # Load saved queries
    st.write("### Load Saved Queries")
    saved_searches = {}
    saved_searches_path = Path("data/saved_searches.json")
    if saved_searches_path.exists():
        try:
            with open(saved_searches_path, "r") as f:
                saved_searches = json.load(f)
        except:
            saved_searches = {}
    
    if saved_searches:
        selected_query = st.selectbox("Select Saved Query", list(saved_searches.keys()))
        if st.button("Load") and selected_query:
            search_manager.current_query = saved_searches[selected_query]
            st.rerun() 