# PubMed Summary App

## Overview

The PubMed Summary App is a Streamlit application designed to help users search PubMed for scientific papers and generate summaries of their findings. It uses the Ollama API with the deepseek-r1:8b model to generate search strategies and summarize search results.

## Features

- **PubMed Integration**: Search PubMed for scientific papers
- **Search Strategy Generation**: Generate PubMed search strategies based on research questions
- **Advanced Search**: Build complex search queries with Boolean operators and field tags
- **Filtering**: Filter search results by year and journal
- **Summarization**: Generate summaries of search results using Ollama
- **Export**: Download summaries as text files and search results as CSV files

## Prerequisites

- [Ollama](https://ollama.ai/) installed locally with the deepseek-r1:8b model
- Email address for PubMed queries

## Installation

1. Clone this repository:
```bash
git clone https://github.com/yourusername/Summary-App.git
cd Summary-App
```

2. Install the required dependencies:
```bash
pip install -r requirements.txt
```

3. Set up environment variables:
```bash
cp .env.example .env
```
Then edit the `.env` file with your PubMed email and other configuration options.

4. Make sure Ollama is running with the deepseek-r1:8b model:
```bash
ollama run deepseek-r1:8b
```

## Usage Guide

1. **Start the application**:
```bash
./start_app.sh
```
   Or alternatively:
```bash
python ssl_bypass.py && streamlit run app.py
```

2. **Define Research Question**:
   Enter your research question and confirm the Ollama model settings.

3. **PubMed Search**:
   - Generate a search strategy based on your research question
   - Or build a search query manually using the advanced search options
   - Enter your PubMed email and execute the search

4. **Filter Results**:
   - Filter the search results by year and journal
   - View the filtered results in a table

5. **Generate Summary**:
   - Click the "Summarize Papers" button to generate a summary of the filtered results
   - The summary will be displayed in the app

6. **Export Results**:
   - Download the summary as a text file
   - Download the search results as a CSV file

## SSL Certificate Issues

The application includes an SSL bypass mechanism to handle certificate verification issues when connecting to PubMed. This is particularly useful on macOS systems where Python may have trouble with SSL certificate verification. The bypass is automatically enabled when you run the application.

## Limitations

Note that while the LLM has been prompted to provide accurate summaries, hallucinations can occur. Always double-check the generated summaries against the original papers. The app is intended as a research aid, not a replacement for critical reading and analysis.

## Acknowledgments

- Built with Streamlit and Ollama
- Uses the Entrez API for PubMed integration
- Inspired by the [Ollama-Review-App](https://github.com/chk-AI/Ollama-Review-App) 