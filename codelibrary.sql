-- 创建分区表
CREATE TABLE genome_data (
    chromosome VARCHAR(2),
    start_pos INT,
    end_pos INT,
    gene_name VARCHAR(50),
    expression_level FLOAT
)
PARTITION BY LIST (chromosome) (
    PARTITION p1 VALUES IN ('1'),
    PARTITION p2 VALUES IN ('2'),
    PARTITION p3 VALUES IN ('3'),
    PARTITION p24 VALUES IN ('X', 'Y')
);


drop table genome_data;
-- 创建分区表
CREATE TABLE genome_data (
    chromosome VARCHAR(10),
    start_pos INT,
    start_posindex INT,
    end_pos INT,
    gene_name VARCHAR(50),
    expression_level FLOAT
)
PARTITION BY LIST COLUMNS(chromosome, start_posindex) (
    PARTITION p0_0 VALUES IN (('chr1',0)),
    PARTITION p0_1 VALUES IN (('chr1',1)),
    PARTITION p1_0 VALUES IN (('chr2',0)),
    PARTITION p1_1 VALUES IN (('chr2',1)),
    PARTITION p2_0 VALUES IN (('chr3',0)),
    PARTITION p2_1 VALUES IN (('chr3',1))
);

-- 在分区表上创建复合索引
CREATE INDEX idx_genome_region
ON genome_data (chromosome, start_pos, end_pos);
