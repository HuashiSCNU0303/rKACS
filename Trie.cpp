#include "Trie.h"
// 应该也没有问题，但是后面换了实现的话应该还要再测试一遍

bool UpdateTrie(Trie* root, Clique& clique, int index)
{
	Trie* p = root;
	bool flag = false;
	for (auto node : clique)
	{
		auto it = p->childs.find(node);
		if (it == p->childs.end())
		{
			flag = true;
			it = p->childs.emplace(node, new Trie()).first;
		}
		// p = p->childs.at(node);
		p = it->second;
		p->index = index;
	}

	return flag; // 标记是否成功插入Trie
}

void GetCandsOnTrie(Trie* root, Clique& results)
{
	if (root->childs.size() == 0)
	{
		return;
	}
	queue<Trie*> Q;
	Q.push(root);
	while (!Q.empty())
	{
		Trie* node = Q.front();
		Q.pop();

		if (node->childs.size() == 0)
		{
			results.push_back(node->index);
		}
		else
		{
			for (auto& data_p : node->childs)
			{
				Q.push(data_p.second);
			}
		}
	}
}

void DestroyTrie(Trie* root)
{
	queue<Trie*> Q;
	Q.push(root);
	while (!Q.empty())
	{
		Trie* node = Q.front();
		Q.pop();

		for (auto& data_p : node->childs)
		{
			Q.push(data_p.second);
		}

		delete node;
	}
}