import pandas as pd
import numpy as np
from matchmaker.dataset import box2d

def test_box2d():
    # Create a sample DataFrame
    df = pd.DataFrame({
        'ra': [1, 2, 3, 4, 5],
        'dec': [10, 20, 30, 40, 50]
    })

    # Define the bounding box coordinates
    box_xxyy = [2, 4, 20, 40]

    # Call the box2d function
    filtered_df = box2d(df, box_xxyy)

    # Define the expected filtered DataFrame
    expected_df = pd.DataFrame({
        'ra': [2, 3, 4],
        'dec': [20, 30, 40]
    })

    # Check if the filtered DataFrame matches the expected DataFrame
    assert filtered_df.equals(expected_df)

# Run the test
test_box2d()import pandas as pd
import numpy as np
from matchmaker.dataset import box2d

def test_box2d():
    # Create a sample DataFrame
    df = pd.DataFrame({
        'ra': [1, 2, 3, 4, 5],
        'dec': [10, 20, 30, 40, 50]
    })
    # Define the bounding box coordinates
    box_xxyy = [2, 4, 20, 40]
    # Call the box2d function
    filtered_df = box2d(df, box_xxyy)
    # Define the expected filtered DataFrame
    expected_df = pd.DataFrame({
        'ra': [2, 3, 4],
        'dec': [20, 30, 40]
    })
    # Check if the filtered DataFrame matches the expected DataFrame
    assert filtered_df.equals(expected_df)

def test_box2d_empty():
    # Create an empty DataFrame
    df = pd.DataFrame(columns=['ra', 'dec'])
    # Define the bounding box coordinates
    box_xxyy = [2, 4, 20, 40]
    # Call the box2d function
    filtered_df = box2d(df, box_xxyy)
    # Check if the filtered DataFrame is empty
    assert filtered_df.empty

def test_box2d_no_filter():
    # Create a sample DataFrame
    df = pd.DataFrame({
        'ra': [1, 2, 3, 4, 5],
        'dec': [10, 20, 30, 40, 50]
    })
    # Define the bounding box coordinates that don't filter any rows
    box_xxyy = [0, 6, 0, 60]
    # Call the box2d function
    filtered_df = box2d(df, box_xxyy)
    # Check if the filtered DataFrame is the same as the original DataFrame
    assert filtered_df.equals(df)

# Run the tests
test_box2d()
test_box2d_empty()
test_box2d_no_filter()