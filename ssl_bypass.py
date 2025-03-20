"""
SSL Bypass Module

This module disables SSL certificate verification globally for the Python session.
It's used to solve SSL certificate verification issues that commonly occur on macOS
and other systems when connecting to external APIs like PubMed.

WARNING: Disabling SSL certificate verification is a security risk and should only
be used for testing or in controlled environments. It makes your application
vulnerable to man-in-the-middle attacks.
"""

import ssl
import urllib.request

# Disable SSL certificate verification globally
ssl._create_default_https_context = ssl._create_unverified_context

# Print confirmation
print("SSL certificate verification has been disabled for this Python session.")
print("This is a security risk and should only be used for testing purposes.")

# Test the connection to PubMed
print("\nTesting connection to PubMed...")
try:
    response = urllib.request.urlopen("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/")
    print(f"Connection successful! Status code: {response.getcode()}")
except Exception as e:
    print(f"Connection failed: {e}")

if __name__ == "__main__":
    print("\nThis script can be run directly before running your app:")
    print("python ssl_bypass.py && streamlit run app.py") 