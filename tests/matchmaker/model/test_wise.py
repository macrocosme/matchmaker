def test_is_AGN():
    # Test case 1: Object is an AGN
    w1w2 = 1.5
    w2w3 = 3.5
    assert is_AGN(w1w2, w2w3) == True

    # Test case 2: Object is not an AGN
    w1w2 = 1.0
    w2w3 = 4.5
    assert is_AGN(w1w2, w2w3) == False

    # Test case 3: Object is on the boundary of AGN classification
    w1w2 = 1.7
    w2w3 = 2.2
    assert is_AGN(w1w2, w2w3) == False

    # Test case 4: Object is outside the AGN color-color box
    w1w2 = 1.5
    w2w3 = 1.5
    assert is_AGN(w1w2, w2w3) == False

# Run the tests
test_is_AGN()