import asyncio
from markdown_it import MarkdownIt
from pyppeteer import launch
import os

# Path to your Markdown file
markdown_file_path = 'docs/output.md'
output_png_path = 'docs/images/toolsflow.png'

# Read the Markdown file
with open(markdown_file_path, 'r', encoding='utf-8') as file:
    markdown_content = file.read()

# Convert Markdown to HTML
md = MarkdownIt()
html_content = md.render(markdown_content)

# Wrap the HTML content in a full HTML document
html_template = f'''
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Markmap</title>
    <style>
        body {{
            margin: 0;
            padding: 0;
            display: flex;
            justify-content: center;
            align-items: center;
            height: 100vh;
        }}
    </style>
</head>
<body>
    {html_content}
</body>
</html>
'''

# Function to export HTML to PNG using Pyppeteer
async def export_html_to_png(html_content, output_path):
    browser = await launch()
    page = await browser.newPage()
    await page.setContent(html_content)
    await page.screenshot({'path': output_path, 'fullPage': True})
    await browser.close()

# Run the export function
asyncio.get_event_loop().run_until_complete(export_html_to_png(html_template, output_png_path))

print(f'Exported {markdown_file_path} to {output_png_path}')
