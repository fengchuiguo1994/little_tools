import mysql.connector
import random
import threading
import time

host = "localhost"
user = "root"
password = "kkltmax8"
spe = "assembly"
precision = 4
header = "header"
nrecord = 5000
binlength = 5000000
database = "test"

# 连接 MySQL 数据库
db = mysql.connector.connect(
    host = host,
    user = user,
    password = password,
    database = database,
    connect_timeout=60
)

cursor = db.cursor()

# 创建表
cursor.execute("""
    CREATE TABLE IF NOT EXISTS your_table5 (
        id INT AUTO_INCREMENT PRIMARY KEY,
        col1 VARCHAR(50),
        col2 INT,
        col3 FLOAT
    )
""")

# 总共插入 1,000,000 行记录
total_rows = 200
batch_size = 5000
num_threads = 5

# 线程函数
def insert_data(start_idx, end_idx):
    try:
        print(2)
        for i in range(start_idx, end_idx):
            # 生成 5,000 行测试数据
            data = [(
                f"value_{j}",
                random.randint(1, 100),
                random.uniform(1.0, 100.0)
            ) for j in range(batch_size)]

            # 在第 500 次循环中故意制造错误
            if i == 499:
                data[0] = (None, None, None)  # 故意制造错误

            # 插入数据
            print(4)
            sql = "INSERT INTO your_table5 (col1, col2, col3) VALUES (%s, %s, %s)"
            print(5)
            cursor.executemany(sql, data)
            print(6)
            db.commit()  # 提交事务
            print(7)
            print(f"Thread {threading.current_thread().name} inserted {batch_size} rows. Total: {(i+1)*batch_size}/{total_rows}")
    except mysql.connector.Error as err:
        print(f"Error occurred: {err}")
        db.rollback()  # 回滚事务

        # 删除创建的表
        cursor.execute("DROP TABLE your_table5")
        db.commit()

        print("Database has been rolled back to the state before creating the table.")
        return

# 创建并启动线程
threads = []
for i in range(num_threads):
    start_idx = i * (total_rows // num_threads)
    end_idx = (i + 1) * (total_rows // num_threads)
    if i == num_threads - 1:
        end_idx = total_rows  # 确保最后一个线程处理完全部剩余数据
    print(start_idx, end_idx)
    t = threading.Thread(target=insert_data, args=(start_idx, end_idx))
    t.start()
    threads.append(t)

# 等待所有线程完成
for t in threads:
    t.join()

db.close()