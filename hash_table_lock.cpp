#include <cstring>
#include <cstdlib>
#include <cstdio>
#include <omp.h>
#include <stdexcept>
#include "hash_table_lock.h"



HashTable::HashTable(size_t size) : table_size(size) {
        table = new Node*[table_size];
        locks = new omp_lock_t[table_size];
        for (size_t i = 0; i < table_size; i++) {
            table[i] = nullptr;
            omp_init_lock(&locks[i]);
        }
    }

HashTable::~HashTable() {
        for (size_t i = 0; i < table_size; i++) {
            omp_destroy_lock(&locks[i]);
            Node* entry = table[i];
            while (entry != nullptr) {
                Node* prev = entry;
                entry = entry->next;
                
               
                CorrectNode* valueNode = prev->values;
                while (valueNode != nullptr) {
                    CorrectNode* valuePrev = valueNode;
                    valueNode = valueNode->next;
                    free(valuePrev->value);  
                    delete valuePrev;        
                }

                free(prev->key);  
                delete prev;      
            }
        }
        delete[] table;
        delete[] locks;
    }

unsigned int HashTable::hashFunction(const wchar_t* key) {
        unsigned int hash = 0;
        while (*key) {
            hash = (hash * 31) + *key;
            key++;
        }
        return hash % table_size;
    }

int HashTable::insert(const wchar_t* key, const wchar_t* value) {
        unsigned int hashIndex = hashFunction(key);
        omp_set_lock(&locks[hashIndex]);
        
        Node* entry = table[hashIndex];
        while (entry != nullptr) {
            if (wcscmp(entry->key, key) == 0) {
               
                if (wcscmp(value, L"") == 0) {
               
                    omp_unset_lock(&locks[hashIndex]);
                    return 1;
                }
                CorrectNode* newValueNode = new CorrectNode;
                newValueNode->value = wcsdup(value);
                
                CorrectNode* current = entry->values;
                CorrectNode* previous = nullptr;
                //insert into the right place the link is order
                while (current != nullptr && wcscmp(current->value, newValueNode->value) < 0) {
                    previous = current;
                    current = current->next;
                }

                
                if (previous == nullptr) {
                    newValueNode->next = entry->values;
                    entry->values = newValueNode;
                } else {
                    
                    newValueNode->next = current;
                    previous->next = newValueNode;
                }

                entry->count++;
                omp_unset_lock(&locks[hashIndex]);
                return 1; 
            }
            entry = entry->next;
        }

        
        Node* newNode = new Node;
        newNode->key = wcsdup(key);
        newNode->count = 0; 

        if (wcscmp(value, L"") != 0) {
            newNode->values = new CorrectNode; 
            newNode->values->value = wcsdup(value); 
            newNode->values->next = nullptr; 
            newNode->count = 1; 
        } else {
            newNode->values = nullptr;
        }
       
        newNode->next = table[hashIndex]; 
        table[hashIndex] = newNode;
        omp_unset_lock(&locks[hashIndex]);
        return 1; 
    }
int HashTable::remove(const wchar_t* key) {
    unsigned int hashIndex = hashFunction(key); 
    omp_set_lock(&locks[hashIndex]);

    Node* entry = table[hashIndex];
    Node* prev = nullptr;

    
    while (entry != nullptr) {
        if (wcscmp(entry->key, key) == 0) {
           
            CorrectNode* valueNode = entry->values;
            while (valueNode != nullptr) {
                CorrectNode* valuePrev = valueNode;
                valueNode = valueNode->next;
                free(valuePrev->value);  
                delete valuePrev;        
            }

            free(entry->key);

            if (prev == nullptr) {
                table[hashIndex] = entry->next;
            } else {
                prev->next = entry->next;
            }

            delete entry; 
            omp_unset_lock(&locks[hashIndex]);
            return 1; 
        }

        prev = entry;
        entry = entry->next;
    }

    
    omp_unset_lock(&locks[hashIndex]);
    return 0; 
}


   
wchar_t** HashTable:: getAllKeys(size_t& count) {
        count = 0; 

        for (size_t i = 0; i < table_size; i++) {
            Node* entry = table[i];
            while (entry != nullptr) {
                count++;
                entry = entry->next;
            }
        }

        
        wchar_t** keys = new wchar_t*[count];
        size_t index = 0; 

        for (size_t i = 0; i < table_size; i++) {
            Node* entry = table[i];
            while (entry != nullptr) {
                keys[index++] = wcsdup(entry->key); 
                entry = entry->next; 
            }
        }

        return keys; 
    }

   
wchar_t** HashTable:: getAllValue(const wchar_t* dictWord, size_t& count) {
        unsigned int hashIndex = hashFunction(dictWord);
        Node* entry = table[hashIndex];

        while (entry != nullptr) {
            if (wcscmp(entry->key, dictWord) == 0) {
                wchar_t** correctWords = new wchar_t*[count]; 

                CorrectNode* valueEntry = entry->values;
                size_t index = 0;
                while (valueEntry != nullptr) {
                    correctWords[index++] = wcsdup(valueEntry->value); 
                    valueEntry = valueEntry->next; 
                }

                return correctWords; 
            }
            entry = entry->next; 
        }

        throw std::runtime_error("Wrong word not found!"); 
    }

int HashTable::getcnt(const wchar_t* dictWord) {
        unsigned int hashIndex = hashFunction(dictWord);
        Node* entry = table[hashIndex];
        if (entry != nullptr && wcscmp(entry->key, dictWord) == 0) {
            return entry->count;
        }
        return -1; 
    }

KeyCountPair* HashTable::getAllKeysAndCounts(size_t& count) {
    count = 0; // Initialize the counter

    // First, count the number of unique keys
    for (size_t i = 0; i < table_size; i++) {
        Node* entry = table[i];
        while (entry != nullptr) {
            count++;
            entry = entry->next;
        }
    }
    // Allocate an array for the key-count pairs
    KeyCountPair* keyCountPairs = new KeyCountPair[count];
    size_t index = 0;

    // Fill the array with keys and their counts
    for (size_t i = 0; i < table_size; i++) {
        Node* entry = table[i];
        while (entry != nullptr) {
            keyCountPairs[index].key = wcsdup(entry->key); // Copy the key
            keyCountPairs[index].count = entry->count;     // Store the count
            index++;
            entry = entry->next;
        }
    }

    return keyCountPairs; // Return the array of key-count pairs
}

void send_serialize(HashTable& hashTable,int rank, int receiver_rank, MPI_Comm comm) {
        MPI_Request request1,request2;
 
        size_t total_size = 0;

        // get the len
        for (size_t i = 0; i < hashTable.table_size; ++i) {
            Node* current = hashTable.table[i];
            while (current) {
                total_size +=((wcslen(current->key) + 1) * sizeof(wchar_t));
                total_size += sizeof(size_t); 
                CorrectNode* correctCurrent = current->values;
                while (correctCurrent) {
                    total_size += ((wcslen(correctCurrent->value) + 1) * sizeof(wchar_t));
                    correctCurrent = correctCurrent->next;
                }
                current = current->next;
            }
        }

        // 发送总大小
        MPI_Isend(&total_size, 1, MPI_UNSIGNED_LONG, receiver_rank, 0, comm, &request1);
        MPI_Wait(&request1, MPI_STATUS_IGNORE); 
        char* buffer = new char[total_size];
        int position=0;

        // create buffer
        for (size_t i = 0; i < hashTable.table_size; ++i) {
            if(hashTable.table[i]==nullptr){
                continue;
            }
            Node* current = hashTable.table[i];
            while (current) {
                size_t keyLen = wcslen(current->key) + 1;

                memcpy(buffer+position, current->key, keyLen*sizeof(wchar_t));
                position+=(keyLen * sizeof(wchar_t));

                
                memcpy(buffer+position, &current->count, sizeof(size_t));
                position += sizeof(size_t) ;

                CorrectNode* correctCurrent = current->values;
                while (correctCurrent) {
                    size_t valueLen = wcslen(correctCurrent->value) + 1;
                    memcpy(buffer+position,  correctCurrent->value, valueLen*sizeof(wchar_t));
                    position+=valueLen * sizeof(wchar_t);
                    correctCurrent = correctCurrent->next;
                }
                current = current->next;
            }
        }

        
        MPI_Isend(buffer, total_size , MPI_BYTE, receiver_rank, 1, comm, &request2);
        MPI_Wait(&request2, MPI_STATUS_IGNORE); 
      
        delete[] buffer;
}

void rec_deserialize(HashTable& hashTable, char* buffer,size_t total_size) {

        size_t position=0;

        while(position<total_size){
            //get the key string
            wchar_t* key_string = reinterpret_cast<wchar_t*>(buffer + position);
            int keylen = wcslen(key_string) + 1; 
            wchar_t* key = new wchar_t[keylen];
            memcpy(key, key_string, keylen * sizeof(wchar_t));
            position += keylen * sizeof(wchar_t);

            //get the cnt of value of the key
            size_t count;
            memcpy(&count, buffer+position, sizeof(size_t));
            position += sizeof(size_t);

            Node* newNode = new Node();
            newNode->key = key;
            newNode->values = nullptr;
            newNode->count = count;
            newNode->next = nullptr;

            for(size_t i=0;i<count;i++){
                wchar_t* current_string = reinterpret_cast<wchar_t*>(buffer + position);
                int len = wcslen(current_string) + 1; 
                wchar_t* value = new wchar_t[len];
                memcpy(value, current_string, len * sizeof(wchar_t)); 
                position += len * sizeof(wchar_t);
                
                CorrectNode* newCorrectNode = new CorrectNode();
                newCorrectNode->value = value;
                newCorrectNode->next = newNode->values;
                newCorrectNode->next = nullptr;

                if (newNode->values == nullptr) {
                    newNode->values = newCorrectNode;
                } else {
                    CorrectNode* last = newNode->values;
                    while (last->next != nullptr) {
                        last = last->next;
                    }
                    last->next = newCorrectNode;
                }
            }

        
            unsigned int index = hashTable.hashFunction(newNode->key);
            omp_set_lock(&hashTable.locks[index]);
            newNode->next = hashTable.table[index];
            hashTable.table[index] = newNode;
            omp_unset_lock(&hashTable.locks[index]);
            
        }
        

        delete[] buffer;
}


