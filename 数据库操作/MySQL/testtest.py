import pymysql
import traceback
from multiprocessing import Pool,Manager,cpu_count
from multiprocessing.managers import BaseManager
import os,sys,time
import random


# 建库建表
def createTable():
    conn=pymysql.connect(
        host='localhost',
        port=3306,
        user='root',
        passwd='kkltmax8',
    )
    cur=conn.cursor()
    sql_database='create database if not exists py3_userinfo default charset utf8 collate utf8_general_ci;'
    sql_table='''
        create table student(
          id int not null auto_increment,
          name varchar(20) not null,
          age int default 0,
          tel varchar(13) unique not null,
          primary key(id)
        )engine=innodb character set utf8;
    '''
    try:
        # 建库
        cur.execute(sql_database)
        conn.select_db('py3_userinfo')
        cur.execute('drop table if exists student;')
        # 建表
        cur.execute(sql_table)
    except Exception as e:
        print(e)
    else:
        print('数据库及表创建成功！')
        cur.close()
        conn.close()


class MyMysql(object):
    def __init__(self):
        # 数据库连接
        self.conn = pymysql.connect(
            host='localhost',
            port=3306,
            user='root',
            passwd='kkltmax8',
            db='py3_userinfo',
            charset='utf8'
        )
        # 游标
        self.cur = self.conn.cursor()

    def executeSql(self,sql):
        try:
            res=self.cur.execute(sql)
            print('执行sql受影响的行数：',res)
            self.conn.commit()
        except:
            traceback.print_exc()

    def quit(self):
        self.cur.close()
        self.conn.commit()
        self.conn.close()


class MyManager(BaseManager):
    pass


def my_Manager():
    m=MyManager()
    m.start()
    return m


# 把myMysql类注册到MyManager管理类中
MyManager.register('MyMysql',MyMysql)


def run(my_sql):
    print('subprocess is ',os.getpid())
    # 造数据
    name = 'ha'+str(random.randint(1,100))+'_'+str(os.getpid())
    age = random.randint(1,100)
    # emil=name+'@qq.com'
    tel='1'+str(random.choice([3,5,7,8]))+str(random.random())[2:11]
    sql="insert into student(name,age,tel) values('%s','%s','%s')"%(name,age,tel)
    my_sql.executeSql(sql)


if __name__ == '__main__':
    createTable()

    manager=my_Manager()
    my_sql=manager.MyMysql()

    print('Parent process ',os.getpid())
    p=Pool(cpu_count())
    n=10
    for i in range(n):
        p.apply(run,args=(my_sql,))
    p.close()
    p.join()
    print('all subprocess done!')
    my_sql.quit()