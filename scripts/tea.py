from os import path
import json
import StringIO as sio
from collections import OrderedDict
import re

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

  _version = "0.1.0"
  _meta = {"invocation":""}

  _tree = "" #TODO internally as an actual tree? convert on write?
  _views = {}

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
  def set_invocation(self, invocation_string):
    self._meta["invocation"] = invocation_string

  def get_invocation(self):
    return self._meta["invocation"]

  def set_version(self, version_string):
    self._version = version_string

  def get_version(self):
    return self._version

  def set_tree(self, tree_string):
    self._tree = tree_string

  def get_tree(self):
    return self._tree

  # ==========================================
  # Modifiers
  # ==========================================

  def add_annotation(self, view_name, edge_id, annotations):
    """ adds an arbitrary number of key-value pairs ("annotations")
        belonging to an edge in the tree ("edge_id")
        to a given named sample ("view_name").
    """
    edge_id = int(edge_id)

    # look up sample, make one if it's not there
    if view_name in self._views.keys():
      if edge_id in self._views[view_name]["annotation"]:
        self._views[ view_name ][ "annotation" ][ edge_id ].update( annotations )
      else:
        self._views[ view_name ][ "annotation" ][ edge_id ] = OrderedDict( sorted(annotations.items()) )
    else:
      self._views[view_name] = {"annotation":
                                  {
                                    edge_id: OrderedDict( sorted(annotations.items()) )
                                  }
                                }

  def annotated_tree( self, view_name, annotation_key, alias_name=None ):
    """ returns an edge-annotated tree, according to the specified annotation_key
    """

    # first build a table of annotations that is more easily accessible
    annotation_lookup = dict()

    if not view_name in self._views.keys():
      raise RuntimeError( "View with name '{}' not found!".format( view_name ) )

    view = self._views[ view_name ]

    for edge_id, annotations in view["annotation"].iteritems():
      if annotation_key in annotations.keys():
        annotation_lookup[ edge_id ] = annotations[ annotation_key ]
      # else:
      #   raise RuntimeError( "Annotation key '{}' not found for edge '{}'!".format( annotation_key, edge_id ) )

    # now parse over the tree string, subsituting {<edge_id>}
    # with [<annotation>] where edge_id has an annotation (erase otherwise)
    tree_array = re.split("[{}]+", self._tree)
    for i in range(1, len(tree_array), 2):
      edge_id = int(tree_array[i])
      if edge_id in annotation_lookup:
        # replace with annotation
        nhx_name = alias_name if alias_name else annotation_key
        tree_array[i] = "[&&NHX:{}={}]".format( nhx_name, annotation_lookup[ edge_id ] )
      else:
        # erase
        tree_array[i] = ""

    return "".join( tree_array )

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
        ("views", o._views),
        ("meta", o._meta),
        ("version", o._version)
      ])
    else:
      return super(TEAJSONEncoder, self).default(o)
