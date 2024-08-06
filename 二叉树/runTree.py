# -*- coding: utf-8 -*-
"""
Created on Fri Jul 31 12:25:00 2020

@author: huang
"""

class Node:
    def __init__(self,right,left,p,color,inter,maxx):
        self.key=inter.low
        self.right=right
        self.left=left
        self.p=p
        self.color=color
        #新增区间属性
        self.inter=inter
        #新增附加信息maxx
        self.maxx=maxx
 
#代表区间的类
class Inter:
    def __init__(self,low,high):
        self.low=low
        self.high=high
 
class tree:
    def __init__(self,root,nil):
        self.root=root
        self.nil=nil
    def tree_insert(self,z):
        y=self.nil
        x=self.root
        while x!=self.nil:
            y=x
            if z.key<x.key:
                x=x.left
            else:
                x=x.right
        z.p=y
        if y==self.nil:
            self.root=z
        elif z.key<y.key:
            y.left=z
        else:
            y.right=z
        z.left=self.nil
        z.right=self.nil
        z.color="RED"
        z.maxx=max(z.inter.high,z.left.maxx,z.right.maxx)
        #红黑树性质维护
        self.rb_insert_fixup(z)
        #更新父结点直到根结点的maxx
        while z.p!=self.nil:
            z.p.maxx=max(z.p.maxx,z.maxx)
            z=z.p
            
    def left_rotate(self,x):
        y=x.right
        x.right=y.left
        if y.left!=self.nil:
            y.left.p=x
        y.p=x.p
        if x.p==self.nil:
            self.root=y
        elif x==x.p.left:
            x.p.left=y
        else:
            x.p.right=y
        y.left=x
        x.p=y
        #左旋导致两个结点的max属性改变，更新如下
        y.maxx=x.maxx
        x.maxx=max(x.left.maxx,x.right.maxx,x.inter.high)
        
 
    def right_rotate(self,y):
        x=y.left
        y.left=x.right
        if x.right!=self.nil:
            x.right.p=y
        x.p=y.p  
        if y.p==self.nil:
            self.root=x
        elif y==y.p.left:
            y.p.left=x
        else:
            y.p.right=x
        x.right=y
        y.p=x
        #右旋导致两个结点的max属性改变，更新如下
        x.maxx=y.maxx
        y.maxx=max(y.right.maxx,y.left.maxx,y.inter.high)
 
    def rb_insert_fixup(self,z):
        while z.p.color=="RED":
            if z.p==z.p.p.left:
                y=z.p.p.right
                if y.color=="RED":
                    z.p.color="BLACK"
                    y.color="BLACK"
                    z.p.p.color="RED"
                    z=z.p.p
                else:
                    if z==z.p.right:
                        z=z.p
                        self.left_rotate(z)
                    z.p.color="BLACK"
                    z.p.p.color="RED"
                    self.right_rotate(z.p.p)
            else:
                y=z.p.p.left
                if y.color=="RED":
                    z.p.color="BLACK"
                    y.color="BLACK"
                    z.p.p.color="RED"
                    z=z.p.p
                else:
                    if z==z.p.left:
                        z=z.p
                        self.right_rotate(z)
                    z.p.color="BLACK"
                    z.p.p.color="RED"
                    self.left_rotate(z.p.p)
        self.root.color="BLACK"
    def inorder_tree_walk(self,x):
        if x!=self.nil:
            self.inorder_tree_walk(x.left)
            print(x.key)
            self.inorder_tree_walk(x.right)
    def tree_search(self,x,k):
        if x==self.nil or k==x.key:
            return x
        if k < x.key:
            return self.tree_search(x.left,k)
        else: return self.tree_search(x.right,k)
 
    def rb_transplant(self,u,v):
        if u.p==self.nil:
            self.root=v
        elif u==u.p.left:
            u.p.left=v
        else:
            u.p.right=v
        v.p=u.p
    def tree_minimum(self,x):
        while x.left!=self.nil:
            x=x.left
        return x
 
    
    def rb_delete(self,z):
        y=z
        y_original_color=y.color
        if z.left==self.nil:
            x=z.right
            self.rb_transplant(z,z.right)
        elif z.right==self.nil:
            x=z.left
            self.rb_transplant(z,z.left)
        else:
            y=self.tree_minimum(z.right)
            y_original_color=y.color
            x=y.right
            if y.p==z:
                x.p=y
            else:
                self.rb_transplant(y,y.right)
                y.right=z.right
                y.right.p=y
            self.rb_transplant(z,y)
            y.left=z.left
            y.left.p=y
            y.color=z.color
            if y_original_color=="BLACK":
                self.rb_delete_fixup(x)       
            
 
    def rb_delete_fixup(self,x):
        while x!=self.root and x.color=="BLACK":
            if x==x.p.left:
                w=x.p.right
                if w.color=="RED":
                    w.color="BLACK"
                    x.p.color="RED"
                    self.left_rotate(x.p)
                    w=x.p.right
                if w.left.color=="BLACK" and w.right.color=="BLACK":
                    w.color="RED"
                    x=x.p
                else:
                    if w.right.color=="BLACK":
                        w.left.color=="BLACK"
                        w.color="RED"
                        self.right_rotate(w)
                        w=x.p.right
                    w.color=x.p.color
                    x.p.color="BLACK"
                    w.right.color="BLACK"
                    self.left_rotate(x.p)
                    x=self.root
            else:
                w=x.p.left
                if w.color=="RED":
                    w.color="BLACK"
                    x.p.color="RED"
                    self.right_rotate(x.p)
                    w=x.p.left
                if w.right.color=="BLACK" and w.left.color=="BLACK":
                    w.color="RED"
                    x=x.p
                else:
                    if w.left.color=="BLACK":
                        w.right.color=="BLACK"
                        w.color="RED"
                        self.left_rotate(w)
                        w=x.p.left
                    w.color=x.p.color
                    x.p.color="BLACK"
                    w.left.color="BLACK"
                    self.right_rotate(x.p)
                    x=self.root
        x.color="BLACK"
 
    def print_tree(self,z):
        if z!=self.nil:
            print(z.key,z.color,"[",z,inter.low,",",z.inter.high,"]",":",end='')
            print("( ",end='')
            print(z.left.key,z.left.color,"   ",end='')
            print(z.right.key,z.right.color,end='')
            print(" )",end='')
    def interval_search(self,i):
        x=self.root
        while x!=self.nil and (i.high<x.inter.low or x.inter.high<i.low):
            if x.left!=self.nil and x.left.maxx>=i.low:
                x=x.left
            else:
                x=x.right
        return x
 
#if __name__=="__main__":
#    inter=Inter(0,0)
#    nil=Node(None,None,None,"BLACK",inter,0)
#    inter=Inter(16,21)
#    root=Node(nil,nil,nil,"BLACK",inter,21)
#    t=tree(root,nil)
#    TT=[8,9,25,30,5,8,15,23,17,19,26,26,0,3,6,10,19,20]
#    for i in range(0,len(TT),2):
#        inter=Inter(TT[i],TT[i+1])
#        z=Node(nil,nil,nil,"RED",inter,0)
#        t.tree_insert(z)
