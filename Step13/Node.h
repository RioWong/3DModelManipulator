#pragma once
#include "MyApp.h"

class Node
{
public:
	Node(void);
	~Node(void);

public:
	Node* doesInclude(int v);
	void append(int v, int f);
	int v;
	int f;
	Node* next;
	Node* tail;

};

