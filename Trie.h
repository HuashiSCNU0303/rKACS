#pragma once

#include "Common.h"

typedef struct Trie
{
	int index;
	flat_hash_map<int, Trie*> childs;
} Trie;

bool UpdateTrie(Trie* root, Clique& clique, int index);

void GetCandsOnTrie(Trie* root, Clique& results);

void DestroyTrie(Trie* root);
