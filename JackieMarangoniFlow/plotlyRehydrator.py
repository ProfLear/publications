import json
import base64
import numpy as np
from bs4 import BeautifulSoup
import plotly.graph_objects as go

def decode_bdata(obj):
    """
    Recursively search for and decode Plotly binary data objects.
    These objects look like: {'dtype': 'f8', 'bdata': '...', 'shape': [10, 10]}
    """
    if isinstance(obj, dict):
        if 'bdata' in obj and 'dtype' in obj:
            # Found a binary data object
            dtype_map = {
                'f8': '<f8',
                'f4': '<f4',
                'i4': '<i4',
                'u4': '<u4',
                'i8': '<i8',
                'u8': '<u8',
            }
            dt = dtype_map.get(obj['dtype'], obj['dtype'])
            raw_data = base64.b64decode(obj['bdata'])
            arr = np.frombuffer(raw_data, dtype=dt)
            
            if 'shape' in obj:
                shape = obj['shape']
                if isinstance(shape, str):
                    shape = [int(x) for x in shape.split(',')]
                elif isinstance(shape, list):
                    shape = [int(x) for x in shape]
                arr = arr.reshape(shape)
            
            return arr.tolist() if isinstance(arr, np.ndarray) else arr
        else:
            return {k: decode_bdata(v) for k, v in obj.items()}
    elif isinstance(obj, list):
        return [decode_bdata(item) for item in obj]
    else:
        return obj

def load_plotly_html(file_path):
    """
    Reads a plotly html file and converts it to a plotly figure object.
    """
    with open(file_path, 'r', encoding='utf-8') as f:
        html_content = f.read()

    soup = BeautifulSoup(html_content, 'html.parser')

    # Find the script tag containing the actual plot data.
    # Usually there are two: one for the Plotly library and one for the data.
    # The data one is typically the last one and contains the specific Plotly.newPlot call.
    target_script = None
    scripts = [tag.string for tag in soup.find_all('script') if tag.string and 'Plotly.newPlot(' in tag.string]
    
    if not scripts:
        raise ValueError("Could not find any script tag containing 'Plotly.newPlot' data.")
    
    # We take the last one as it's most likely the one with the data call,
    # and not the one in the library source code.
    target_script = scripts[-1]
    
    # Locate the "Plotly.newPlot(" signature execution point
    start_idx = target_script.find('Plotly.newPlot(')
    if start_idx == -1:
         raise ValueError("Could not locate 'Plotly.newPlot(' in script.")

    # The arguments are: (elementId, dataArray, layoutObject, configObject)
    # We skip the elementId (first argument) by finding the first comma.
    first_comma = target_script.find(',', start_idx)
    if first_comma == -1:
        raise ValueError("Could not find arguments in Plotly.newPlot call.")
        
    # Find the start of the data array
    data_start = target_script.find('[', first_comma)
    if data_start == -1:
        raise ValueError("Could not locate the start of Plotly data array.")
    
    # Decode the data list using raw_decode to handle the following arguments
    decoder = json.JSONDecoder()
    data, data_end = decoder.raw_decode(target_script[data_start:])
    
    # Recursively decode binary data if present
    data = decode_bdata(data)

    # Look for the next opening brace '{' after the data block (marks the layout object)
    layout_start = target_script.find('{', data_start + data_end)
    if layout_start == -1:
        # Layout might be missing or default, but Plotly.newPlot usually has it
        layout = {}
    else:
        # Decode the layout dictionary
        layout, _ = decoder.raw_decode(target_script[layout_start:])
        layout = decode_bdata(layout)

    # Reconstitute the figure object
    return go.Figure(data=data, layout=layout)

if __name__ == "__main__":
    import sys
    import os

    # Example usage
    example_file = "BarriersBLBottom.html"
    if len(sys.argv) > 1:
        example_file = sys.argv[1]
    
    if not os.path.exists(example_file):
        print(f"File not found: {example_file}")
        sys.exit(1)

    try:
        fig = load_plotly_html(example_file)
        print(f"Successfully rehydrated figure from {example_file}")
        print(f"Number of traces: {len(fig.data)}")
        
        # Print some info about the first trace
        if len(fig.data) > 0:
            trace = fig.data[0]
            print(f"First trace type: {trace.type}")
            if hasattr(trace, 'z') and trace.z is not None:
                z = np.array(trace.z)
                print(f"First trace 'z' shape: {z.shape}")
        
        # If you wanted to show it:
        # fig.show()
        
    except Exception as e:
        print(f"Error: {e}")
        import traceback
        traceback.print_exc()
