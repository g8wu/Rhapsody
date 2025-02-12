import pandas as pd
import os
from bs4 import BeautifulSoup
from io import StringIO

# Input folder, output CSV file
directory = "C:/Users/gio8w/OneDrive - University of California, San Diego Health/Shared Documents - Rhapsody/Pipeline Summaries"
output = "C:/Users/gio8w/OneDrive - University of California, San Diego Health/Shared Documents - Rhapsody/Pipeline Summaries/OUTPUT.csv"

with open(output, 'w', newline='', encoding='utf-8') as csvfile:
    csvwriter = csv.writer(csvfile)
    
    # Loop through all files in the directory
    for filename in os.listdir(directory):
        if filename.endswith('.html'):
            filepath = os.path.join(directory, filename)
            
            # Write file title
            csvwriter.writerow([filename])
            
            # Load HTML content
            with open(filepath, 'r', encoding='utf-8') as file:
                html_content = file.read()
            
            # Parse HTML content
            soup = BeautifulSoup(html_content, 'html.parser')
            
            # Find all tables
            tables = soup.find_all('table')  # Use 'table' tag to find all tables
            if not tables:
                print(f"No tables found in {filename}")
                continue
            
            table_data = []
            for table in tables:
                # Read the table into a DataFrame
                df = pd.read_html(StringIO(str(table)), flavor='html5lib')[0]
                # Append the DataFrame to the list
                table_data.append(df)
            
            # Concatenate all DataFrames into a single DataFrame
            combined_df = pd.concat(table_data, ignore_index=True)
            
            # Write the combined DataFrame to the CSV file
            combined_df.to_csv(csvfile, index=False)
            csvwriter.writerow([])

