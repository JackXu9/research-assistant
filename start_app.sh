#!/bin/bash

# Check if .env file exists
if [ ! -f .env ]; then
    echo "Creating .env file from .env.example..."
    cp .env.example .env
    echo "Please edit .env file with your Groq API key and PubMed email"
    exit 1
fi

# Check if GROQ_API_KEY is set
if ! grep -q "GROQ_API_KEY=" .env; then
    echo "Please set your GROQ_API_KEY in the .env file"
    exit 1
fi

# Check if PUBMED_EMAIL is set
if ! grep -q "PUBMED_EMAIL=" .env; then
    echo "Please set your PUBMED_EMAIL in the .env file"
    exit 1
fi

# Run the app
streamlit run app.py 