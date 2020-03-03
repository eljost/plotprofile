def to_list(inp):
    if isinstance(inp, str):
        inp = [inp, ]
    elif not isinstance(inp, list):
        inp = list(inp)
    return inp
