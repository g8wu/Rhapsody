import os
import csv
from bs4 import BeautifulSoup
from openpyxl import Workbook

#directory = "C:/Users/gio8w/OneDrive - University of California, San Diego Health/Shared Documents - Rhapsody/Pipeline Summaries"
# Home pc
directory = "C:/Users/Ginny/OneDrive - University of California, San Diego Health/Shared Documents - Rhapsody/Pipeline Summaries"
output = directory + "/1_Pipeline_Summaries.csv"
tables_data = []

# Loop through all files in the directory
for filename in os.listdir(directory):
    if filename.endswith('.html'):
        print("Processing file: " + filename)
        filepath = os.path.join(directory, filename)
        
        # Load and parse HTML text
        with open(filepath, 'r', encoding='utf-8') as file:
            html_content = file.read()
        soup = BeautifulSoup(html_content, 'lxml')
        
        # Find all tables
        tables = soup.find_all('table')
        file_tables_data = []
        for table in tables:
            rows = table.find_all('tr')
            
            # Collect headers and data rows
            table_data = []
            headers = [th.text.strip() for th in rows[0].find_all('th')]
            if headers:
                table_data.append([filename])    # Add filename
                table_data.append(headers)
            for row in rows[1:]:
                data = [td.text.strip() for td in row.find_all('td')]
                table_data.append(data)

            file_tables_data.append(table_data)
        tables_data.append(file_tables_data)

# Flatten and transpose
transposed_data = []
for file_tables in tables_data:
    max_rows = max(len(table) for table in file_tables)
    max_columns = max(len(row) for table in file_tables for row in table)
    for row_index in range(max_rows):
        new_row = []
        for table in file_tables:
            if row_index < len(table):
                row = table[row_index]
                new_row.extend(row)
                if len(row) < max_columns:
                    new_row.extend([''] * (max_columns - len(row)))
            else:
                new_row.extend([''] * max_columns)
        transposed_data.append(new_row)

# Write the transposed data to the CSV file
with open(output, 'w', newline='', encoding='utf-8') as csvfile:
    csvwriter = csv.writer(csvfile)
    for row in transposed_data:
        csvwriter.writerow(row)

print("Scrape successful!")
