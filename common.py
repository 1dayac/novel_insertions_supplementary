
class Insertion(object):
    def __init__(self, node, pos1, pos2, rc, ref_id, ref_pos, anchor1, anchor2):
        self.node = node
        self.pos1 = pos1
        self.pos2 = pos2
        self.rc = rc
        self.ref_id = ref_id
        self.ref_pos = ref_pos
        self.anchor1 = anchor1
        self.anchor2 = anchor2

    def __eq__(self, other):
        return abs(self.ref_pos - other.ref_pos) < 100

    def size(self):
        return abs(self.pos2 - self.pos1) - 1


