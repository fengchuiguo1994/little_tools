#include <stdio.h>

int main() {
    FILE *file;
    unsigned char buffer[100];
    size_t bytesRead;

    // 以二进制读模式打开文件
    file = fopen("test1.e500.cpu.cluster", "rb");

    if (file == NULL) {
        printf("无法打开文件\n");
        return 1;
    }

    // 逐块读取并以十六进制方式打印文件内容
    while ((bytesRead = fread(buffer, 1, sizeof(buffer), file)) > 0) {
        for (size_t i = 0; i < bytesRead; i++) {
            printf("%s", buffer[i]); // 以十六进制方式打印每个字节的数值
        }
    }

    // 关闭文件
    fclose(file);

    return 0;
}