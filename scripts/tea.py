from os import path
import json
import StringIO as sio
from collections import OrderedDict

try:
    import textwrap
    textwrap.indent
except AttributeError:  # undefined function (wasn't added until Python 3.3)
    def indent(text, amount, ch=' '):
        padding = amount * ch
        return ''.join(padding+line for line in text.splitlines(True))
else:
    def indent(text, amount, ch=' '):
        return textwrap.indent(text, amount * ch)


class TEA:
  """Class holding Tree Edge Anntations. Ensures adherence to the related file format specifications"""

  # ==========================================
  # Member Variables
  # ==========================================

  _version = "0.0.1"
  _meta = {"invocation":""}

  _tree = "" #TODO internally as an actual tree? convert on write?
  _samples = []

  # ==========================================
  # Constructor
  # ==========================================
  def __init__(self, tree=None, invocation=None, version=None):
    if tree != None:
      self._tree = tree

    if invocation != None:
      self._meta["invocation"] = invocation

    if version != None:
      self._version = version

  # ==========================================
  # Getter/Setter
  # ==========================================
  def invocation(self, invocation_string):
    self._meta["invocation"] = invocation_string

  def invocation(self):
    return self._meta["invocation"]

  def version(self, version_string):
    self._version = version_string

  def version(self):
    return self._version

  def tree(self, tree_string):
    self._tree = tree_string

  def tree(self):
    return self._tree

  # def sample(self, name):
  #   # find sample with specified name

  #   # return it

  # ==========================================
  # Modifiers
  # ==========================================

  def add_annotation(self, sample_name, edge_id, annotations):
    """ adds an arbitrary number of key-value pairs ("annotations")
        belonging to an edge in the tree ("edge_id")
        to a given named sample ("sample_name").
    """

    # look up sample, make one if it's not there
    sample_found = False

    for s in self._samples:
      if s["name"] == sample_name:
        sample_found = True
        edge_found = False
        for a in s["annotation"]:
          if a["edge"] == edge_id:
            edge_found = True
            a.update( annotations )
        if not edge_found:
          # annotations["edge"] = edge_id
          s["annotation"].append(
            OrderedDict( [("edge", edge_id)] + sorted(annotations.items()) )
          )

    if not sample_found:
      annotations["edge"] = edge_id
      self._samples.append( { "name": sample_name,
                              "annotation": [
                                OrderedDict( [("edge", edge_id)] + sorted(annotations.items()) )
                              ]})

  # ==========================================
  # Output
  # ==========================================

  def to_file(self, file_path):
    with open(file_path, "w+") as f:
      f.write( json.dumps( self, cls=TEAJSONEncoder, indent=2 ) )

  def to_stream(self, stream):
    stream.write( json.dumps( self, cls=TEAJSONEncoder, indent=2 ) )

class TEAJSONEncoder(json.JSONEncoder):
  def default(self, o):
    if isinstance(o, TEA):
      return OrderedDict([
        ("tree", o._tree),
        ("samples", o._samples),
        ("meta", o._meta),
        ("version", o._version)
      ])
    else:
      return super(TEAJSONEncoder, self).default(o)
