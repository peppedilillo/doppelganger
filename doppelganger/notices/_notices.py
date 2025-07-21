from pathlib import Path
import xml.etree.ElementTree as Node


def _safe_put(k, v, d: dict):
    """
    Puts value `v` into dictionary `d` with key `k`.
    If key `k` exists, `d[k]` is transformed into a list and `v` is appended to it.
    """
    if k in d:
        if not isinstance(d[k], list):
            d[k] = [d[k]]
        d[k].append(v)
    else:
        d[k] = v


def _parse_tree(node: Node) -> dict:
    """
    DFS implementation of `parse_notice`.
    Parses a single ElementTree node into a dictionary.
    """

    def visit(n: Node):
        """Helper. Node operations, returns a key-value pair."""
        k, v = n.tag, {}
        v.update(**{f"@{k_}": v_ for k_, v_ in n.attrib.items()})
        if n.text is None:
            return k, v
        if not n.text.strip():
            return k, v
        v["#text"] = n.text
        return k, v

    def dfs(n: Node, o: dict):
        """Helper. Depth-first traversal from node `n`, inserts into `o`."""
        if n is None:
            return
        k, v = visit(n)
        _safe_put(k, v, o)
        for child in n:
            dfs(child, v)

    out = {}
    dfs(node, out)
    return out


def parse_notice(content: str) -> dict:
    """
    Parses a VOE notice string into a nested dictionary.
    XML tags are keys to this dictionary.
    Attributes are also keys, but they are noted with '@' prefix.
    Text within tags is stored with key '#text'.
    """
    root = Node.XML(content)
    out = {}
    # we peel out the root node which brings no useful information here
    out.update(_parse_tree(root))
    assert len(out) == 1
    return out.pop(*out.keys())


def parse_notice_file(path: str | Path) -> dict:
    """
    Parses a VOE notice file into a nested dictionary.
    XML tags are keys to this dictionary.
    Attributes are also keys, but they are noted with '@' prefix.
    Text within tags is stored with key '#text'.
    """
    with open(path, "r") as f:
        return parse_notice(f.read())


def _get_nested_key(
    notice: dict,
    name: str,
    key: str,
) -> dict:
    """
    Returns first nested 'key' matching `name` in notice.
    Returns empty dict if no match is found or if notice does not contain `key`.
    """
    if notice == {} or key not in notice:
        return {}
    param_list = notice[key]
    if not isinstance(param_list, list):
        param_list = [param_list]
    for param in param_list:
        if param["@name"] == name:
            return param
    return {}


def get_param(notice: dict, name: str) -> dict:
    """
    Returns first parameter matching `name` in notice.
    Returns empty dict otherwise.

    Usage: > get_param(get_group(notice["What"], "tag-name"), s).get("@value", None))
    """
    return _get_nested_key(notice, name, "Param")


def get_group(notice: dict, name: str) -> dict:
    """
    Returns first group matching `name` in notice.
    Returns empty dict otherwise.

    Usage: > get_param(get_group(notice["What"], "tag-name"), s).get("@value", None))
    """
    return _get_nested_key(notice, name, "Group")
