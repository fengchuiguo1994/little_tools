class Node:
    """
    二叉树节点类
    """
    def __init__(self,right,left,p,color,inter,maxx):
        self.key = inter.low # 区间的起点
        self.right = right # 右节点
        self.left = left # 左节点
        self.p = p
        self.color = color # 红黑标记
        self.inter = inter # 区间类
        self.maxx = maxx # 以当前节点为根的子树中所有区间的端点的最大值
    def __str__(self):
        return self.key
    def __repr__(self):
        return self.__str__()
class Inter:
    """
    表示区间的类
    """
    def __init__(self,low,high):
        self.low = low # 区间起点
        self.high = high # 区间终点
class tree:
    """
    二叉红黑树
    """
    def __init__(self,root,nil):
        self.root = root
        self.nil = nil
    def tree_insert(self,z):
        y = self.nil
        x = self.root
        while x != self.nil:
            y = x
            if z.key < x.key:
                x = x.left
            else:
                x = x.right
        z.p = y
        if y == self.nil:
            self.root = z
        elif z.key < y.key:
            y.left = z
        else:
            y.right = z
        z.left = self.nil
        z.right = self.nil
        z.color = "RED"
        z.maxx = max(z.inter.high,z.left.maxx,z.right.maxx)
        self.rb_insert_fixup(z)
        while z.p != self.nil:
            z.p.maxx = max(z.p.maxx,z.maxx)
            z = z.p