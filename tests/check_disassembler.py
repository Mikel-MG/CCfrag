import sys

sys.path.append("..")

from CCfrag import Divider

list_parameters = [
    ["default", "ok", [20, 10]],
    ["default", "ok", [15, 10]],
    ["default", "ok", [10, 5]],
    ["default", "ok", [10, 0]],
    ["big L limit (no error expected)", "ok", [49, 10]],
    ["big L exception", "error", [50, 0]],
    ["small L exception", "error", [7, 5]],
    ["big O exception", "error", [23, 23]],
]

list_results = []

for concept, expected_result, [L, O] in list_parameters:
    print("#########################################")
    result = None
    print(concept)
    divider = Divider(L, O, 2, "")

    try:
        divider.debug_check()
        result = "ok"
    except Exception as e:
        print("Error was found")
        print(e)
        result = "error"

    test_status = expected_result == result
    list_results.append(test_status)

    print("#########################################")

print(f"{sum(list_results)} out of {len(list_results)} tests were succesful")
