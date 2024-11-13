#ifndef HASHTABLE_H
#define HASHTABLE_H

#include <cstring>
#include <cstdlib>
#include <cstdio>
#include <omp.h>
#include <stdexcept>
#include <mpi.h>

struct KeyCountPair {
    wchar_t* key;    // Pointer to the key (wrong word)
    size_t count;    // Count of correct words associated with the key
};

    // 错误单词节点结构
    struct CorrectNode {
        wchar_t* value;                  // 正确的单词
        CorrectNode* next;              // 指向下一个正确单词节点
    };

    // 错别字节点结构
    struct Node {
        wchar_t* key;                   // 错别字
        CorrectNode* values;           // 指向正确单词链表的指针
        Node* next;                  // 指向下一个错别字节点
        size_t count;                // 记录正确单词的数量
    };

// 哈希表结构
struct HashTable {
    Node** table;                // 哈希表数组
    size_t table_size;           // 哈希表大小
    omp_lock_t* locks;           // 每个槽位的锁数组

    // 构造函数
    HashTable(size_t size);
    
    // 析构函数
    ~HashTable();

    // 哈希函数
    unsigned int hashFunction(const wchar_t* key);

    // 插入错别字和对应的正确单词
    int insert(const wchar_t* key, const wchar_t* value);

    // 获取所有错别字的 key
    wchar_t** getAllKeys(size_t& count);

    // 获取指定错别字的所有正确单词
    wchar_t** getAllValue(const wchar_t* dictWord, size_t& count);
    int remove(const wchar_t* key);

    // 获取指定错别字的正确单词数量
    int getcnt(const wchar_t* dictWord);
    KeyCountPair* getAllKeysAndCounts(size_t& count);
};
void send_serialize(HashTable& hashTable,int rank, int receiver_rank, MPI_Comm comm);
void rec_deserialize(HashTable& hashTable, char* buffer,size_t total_size);
#endif // HASHTABLE_H
