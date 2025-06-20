import mysql.connector
import random

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
    database = database
)

cursor = db.cursor()

try:
    # 创建表
    cursor.execute("""
        CREATE TABLE IF NOT EXISTS your_table1 (
            id INT AUTO_INCREMENT PRIMARY KEY,
            col1 VARCHAR(50) NOT NULL,
            col2 INT NOT NULL,
            col3 FLOAT NOT NULL
        )
    """)

    # 总共插入 1,000,000 行记录
    total_rows = 1000000
    batch_size = 5000

    for i in range(int(total_rows / batch_size)):
        # 生成 5,000 行测试数据
        data = [(
            f"value_{j}",
            random.randint(1, 100),
            random.uniform(1.0, 100.0)
        ) for j in range(batch_size)]

        # 在第 500 次循环中故意制造错误

        # 插入数据
        sql = "INSERT INTO your_table1 (col1, col2, col3) VALUES (%s, %s, %s)"
        cursor.executemany(sql, data)
        db.commit()  # 提交事务

        print(f"Inserted {batch_size} rows. Total: {(i+1)*batch_size}/{total_rows}")

except mysql.connector.Error as err:
    print(f"Error occurred: {err}")
    db.rollback()  # 回滚事务

    # 删除创建的表
    cursor.execute("DROP TABLE your_table1")
    db.commit()

    print("Database has been rolled back to the state before creating the table.")

finally:
    db.close()


logger.info("{0} database created fail".format(args.assembly))