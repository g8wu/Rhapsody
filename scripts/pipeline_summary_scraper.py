import os
import csv
from bs4 import BeautifulSoup

# Directory containing HTML files
directory = "C:/Users/gio8w/OneDrive - University of California, San Diego Health/Shared Documents - Rhapsody/Pipeline Summaries"
output = "C:/Users/gio8w/OneDrive - University of California, San Diego Health/Shared Documents - Rhapsody/Pipeline Summaries/OUTPUT.csv"

# Open a CSV file to write the data
with open(output, 'w', newline='', encoding='utf-8') as csvfile:
    csvwriter = csv.writer(csvfile)
    
    # Write the header row (if any)
    header_written = False
    
    # Loop through all files in the directory
    for filename in os.listdir(directory):
        if filename.endswith('.html'):
            filepath = os.path.join(directory, filename)
            
            # Load the HTML content
            with open(filepath, 'r', encoding='utf-8') as file:
                html_content = file.read()
            
            # Parse the HTML content
            soup = BeautifulSoup(html_content, 'html.parser')
            
            # Find the data you want to scrape
            table = soup.find('table')
            if table:
                rows = table.find_all('tr')
                
                # Write the header row (if not written yet)
                if not header_written and rows:
                    header = [th.text.strip() for th in rows[0].find_all('th')]
                    csvwriter.writerow(header)
                    header_written = True
                
                # Write the data rows
                for row in rows[1:]:
                    data = [td.text.strip() for td in row.find_all('td')]
                    csvwriter.writerow(data)
