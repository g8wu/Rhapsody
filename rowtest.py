import os
from bs4 import BeautifulSoup
import csv
directory = "C:/Users/Ginny/OneDrive - University of California, San Diego Health/Shared Documents - Rhapsody/Pipeline Summaries"
output = directory + "/test.csv"
html_content = '''
<table class="table" style="margin-bottom: 0;">
<tr>
<th class="parentCell border-top-0">
<div data-placement="top" data-toggle="tooltip" title="Type of bioproduct (mRNA, AbSeq, ATAC)">Bioproduct Type</div>
</th>
<th class="parentCell border-top-0">
<div data-placement="top" data-toggle="tooltip" title="Number of filtered read pairs aligned to bioproduct type">Aligned Reads By Type</div>
</th>
<th class="parentCell border-top-0">
<div data-placement="top" data-toggle="tooltip" title="Average number of reads representing the molecules detected in each cell">Mean Reads per Cell</div>
</th>
<th class="parentCell border-top-0">
<div data-placement="top" data-toggle="tooltip" title="Average number of molecules detected per cell label">Mean Molecules per Cell</div>
</th>
</tr>
<tr>
<td class="parentCell" style="padding-right:3.125em">mRNA </td>
<td class="parentCell" style="padding-right:3.125em">1,378,177,625 </td>
<td class="parentCell" style="padding-right:3.125em">32,451.23 </td>
<td class="parentCell" style="padding-right:3.125em">2,290.4 </td>
<tr>
<tr>
<td class="parentCell" style="padding-right:3.125em">AbSeq </td>
<td class="parentCell" style="padding-right:3.125em">836,057 </td>
<td class="parentCell" style="padding-right:3.125em">4.23 </td>
<td class="parentCell" style="padding-right:3.125em">3.79 </td>
<tr>
<tr>
<td class="parentCell" style="padding-right:3.125em">Sample Tags </td>
<td class="parentCell" style="padding-right:3.125em">44,730,780 </td>
<td class="parentCell" style="padding-right:3.125em">791.2 </td>
<td class="parentCell" style="padding-right:3.125em">- </td>
<tr>
</tr></tr></tr></tr></tr></tr></table>
'''

soup = BeautifulSoup(html_content, 'lxml')


with open(output, 'w', newline='', encoding='utf-8') as csvfile:
    csvwriter = csv.writer(csvfile)
    
    # Find the table
    table = soup.find('table', {'class': 'table'})
    
    # Extract headers
    headers = [th.get_text(strip=True) for th in table.find_all('th')]
    csvwriter.writerow(headers)
    
    # Extract rows
    rows = table.find_all('tr')[1:]  # Skip the header row
    for row in rows:
        data = [td.get_text(strip=True) for td in row.find_all('td')]
        if data:  # Only write rows that have data
            csvwriter.writerow(data)

print("------------------Scraping done!------------------")