import asyncio
from pyppeteer import launch
import os

async def convert_html_to_png(html_file, output_file):
    browser = await launch()
    page = await browser.newPage()

    # Convert the HTML file path to an absolute path
    html_file_path = os.path.abspath(html_file)
    # Ensure the path is correctly formatted for the file URL
    file_url = f'file://{html_file_path}'

    await page.goto(file_url, waitUntil='networkidle0')
    element = await page.querySelector('#container')
    bounding_box = await element.boundingBox() # type: ignore
    await page.screenshot({
        'path': output_file,
        'clip': {
            'x': bounding_box['x'], # type: ignore
            'y': bounding_box['y'], # type: ignore
            'width': bounding_box['width'], # type: ignore
            'height': bounding_box['height'] # type: ignore
        }
    })
    await browser.close()

html_file = 'docs/toolsflow.html'  # Replace with the path to your HTML file
output_file = 'output/markmap.png'  # Replace with the desired output PNG path

# Create the output directory if it doesn't exist
os.makedirs(os.path.dirname(output_file), exist_ok=True)

asyncio.get_event_loop().run_until_complete(convert_html_to_png(html_file, output_file))
