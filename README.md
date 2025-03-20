# Your Research Assistant

## Overview

Your Research Assistant is a powerful Streamlit application designed to help researchers efficiently search PubMed for scientific papers and generate intelligent summaries of their findings. It leverages the QwQ-32B reasoning model through the Groq Cloud API to provide high-quality search strategies, comprehensive summaries, and interactive literature analysis.

## Features

- **PubMed Integration**: Advanced search functionality for scientific papers
- **AI-Powered Search Strategy**: Generate optimized PubMed search strategies based on your research questions
- **Advanced Search Options**: 
  - Multiple article types selection
  - Multi-language support
  - Date range filtering
  - Journal filtering
- **Smart Filtering**: Filter search results by year and journal
- **AI Analysis**:
  - Generate comprehensive summaries using QwQ-32B reasoning model
  - "Ask the Literature" feature for specific questions about your papers
  - Support for both full dataset and latest 10 papers analysis
- **Export Options**: 
  - Download summaries as text files
  - Export search results as CSV files
  - Save literature analysis results

## Prerequisites

- Groq Cloud API key (obtain from [console.groq.com/keys](https://console.groq.com/keys))
- Email address for PubMed queries (required by NCBI)

## Installation

1. Clone this repository:
```bash
git clone https://github.com/JackXu9/research-assistant.git
cd research-assistant
```

2. Install the required dependencies:
```bash
pip install -r requirements.txt
```

3. Set up environment variables:
```bash
cp .env.example .env
```
Then edit the `.env` file with your Groq API key and PubMed email.

## Usage Guide

1. **Start the application**:
```bash
streamlit run app.py
```

2. **Configure API Access**:
   - Enter your Groq API key
   - Provide your email for PubMed searches

3. **Research Process**:
   a. **Define Research Question**:
      - Enter your research question
      - Get AI-generated search strategy
   
   b. **PubMed Search**:
      - Use the generated strategy or create your own query
      - Apply filters (article types, languages, date range)
      - Execute the search
   
   c. **Review Results**:
      - View papers grouped by year
      - Expand/collapse paper details
      - Toggle abstracts for detailed review
   
   d. **Generate Insights**:
      - Create full summaries or analyze latest 10 papers
      - Ask specific questions about the literature
      - Download generated content

## Features in Detail

### 1. Search Capabilities
- Multiple article types (Clinical Trials, Reviews, Meta-Analyses, etc.)
- Support for 10 languages including English, French, German, Spanish, etc.
- Custom date range selection
- Advanced query building with automatic filter integration

### 2. AI Analysis
- **Summary Generation**: Comprehensive analysis of selected papers
- **Latest Papers Focus**: Option to analyze most recent 10 papers
- **Interactive Q&A**: Ask specific questions about the selected papers
- **Context-Aware**: Summaries consider your research question

### 3. User Interface
- Clean, intuitive layout
- Year-based paper organization
- Expandable paper details
- Easy-to-use filtering system

## Limitations

- Free tier allows analysis of up to 10 papers per batch
- Additional papers require Groq Cloud subscription
- As with all AI systems, outputs should be verified against source materials

## Security Note

The application handles API keys securely and does not store them permanently. Each user must provide their own Groq API key.

## Support

For issues or questions, please contact: jack.junchi.xu@regionh.dk

## Acknowledgments

- Built with Streamlit
- Powered by Groq Cloud and QwQ-32B reasoning model
- Uses the Entrez API for PubMed integration 