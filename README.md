# 3DModelManipulator

A 3D model manipulator based on ridge-valley lines, and realized with the help of OpenSceneGraph.

## 使用说明

1. 双击3DModelManipulator.exe即可打开程序；
2. 按F键退出全屏，OSG窗口为操作窗口，CMD窗口为状态显示窗口；
3. 在操作窗口鼠标左键点击，选择要操作的点；
4. 程序会找到距离点击点最近的谷脊点，并生成拖拽器；
5. 用鼠标通过拖拽器进行交互式操作；
6. 拖拽器移动到合适位置时，按R键获得移动向量，并更新模型显示；
7. 上述步骤（3-6）可重复进行；
8. 关闭操作窗口，程序将自动保存模型为OBJ和OSG文件；
9. 保存文件完毕后，关闭CMD窗口即可。

## 实现步骤

Step1
实现点的交互式拖动

fixedSizeInScreen的作用，保存文件看一下结点层次

中间过程使用double，更准确

使用点或几何体的区别

---

Step2
实现用鼠标点选模型中的点

	pick()
	handle()
	apply()

遍历所有的点

vector，迭代器的使用？

求距离的公式存在错误，还没更改！

---

Step3
实现选中的点可以交互式操作

在选中的点上生成拖拽器

拖拽器尺寸太大，最小距离的点不准

添加调整拖拽器大小功能。没必要

---

Step4
实现点的移动向量的获取

从哪个函数里获取呢？
在run的过程中有矩阵向量

拖拽器在操控指定模型的时候，通过DraggerCallback来更新实现变换。

traverse递归函数是场景中的每个节点都会调用的，因此在这个函数里，当是事件遍历的时候调
用了内置的handle方法。实现拖拽动作。
这个思路又给我们一些启示，与场景的交互，不一定要继承EventHandler，通过traverse同样可
以实现。这些需要深入理解在事件遍历的时候osg都干了些什么。
addTransformUpdating这个方法，通过回调实现变换。

考虑通过某个按键激活，从而获得移动向量。

到最里层才可以获得移动向量！

按键Refresh，获得向量。

---

Step5
实现给定向量点会即时移动

移动可以通过删除新建完成。

可以移动了，只是嵌套层次有点小问题

---

Step6
实现拖动一个点周围的点的位置也随着变化

先读取模型，再实时更新模型（暂存）

---

Step7
读取PLY2文件

	Group {
		Geode {
			Geometry {
				Primitivesets 1 {
					DrawElementsUInt TRIANGLES 115332 {}
					VertexArray Vec3Array 19611 {}
					NormalBinding PER_VERTEX
					NormalArray Vec3Array 19611 {}
				}
			}
		}
	}

osg读ply

读ply在vertexData.cpp的244行

读中进行的操作

	osg::Geometry* geom  =  new osg::Geometry;
	geom->setVertexArray(_vertices.get());
	geom->setNormalArray(_normals.get());
	geom->setNormalBinding(osg::Geometry::BIND_PER_VERTEX);
	geom->addPrimitiveSet(_triangles.get());  //行448
	// set flage true to activate the vertex buffer object of drawable
	geom->setUseVertexBufferObjects(true);
	osg::Geode* geode = new osg::Geode;
	geode->addDrawable(geom);
	return geode;

11535

fscanf()读浮点数的问题

osg::Vec3Array的声明要注意

光线还是有问题！

终于可以完美的实现导入ply2啦。

---

Step8
模型的实时更新回调，

改变部分点的坐标，能够即时显示。

在操作时，删除原结点，新建几何体。

Group, Node, Drawable, Geometry, get

思路明白，暂存，继续Step6

---

Step9

没有功能。
可能是想只能判断画布数目。

---

Step10

只有修改点的大小。
可能是想直接控制像素点，但是失败了。

---

Step11
单文档

不明白为什么会闪烁
单文档不新建窗口问题
单窗口失败，创建多文档窗口。

---

Step12
多文档

SendStrToOutputWnd还有问题，可能是类型的问题！
要放到该放的地方。

将用OpenGL写的提取谷脊线的部分移植到OSG
将cmd移植到MFC

交互式操作无法进行，总是意外停止。
注释掉也不可以？
mViewer->addEventHandler(new PickHandler(mRoot.get()));的位置有问题！
已经完美的实现了！

---

Step13
一个失败的结合体。

提取谷脊线主要集中在Doc类，
交互式变形主要集中在View类。

单文件版。

位于void MeshData::computePrincipal()

	void MeshData::getConsistent1Ring(int index, int *&nei, int &nei_N)
	{
		nei_N = degree_v[index];
		nei = vertex_link_v[index];
	}

void MeshData::generateVertexLink()得到degree_v
果然在这。

	sortVertexLink(link[i], link_N[i], vertex_link_v[i], vertex_link_f[i], degree_v[i], degree_f[i]);

## TODO

1. MESH三角形面片的直接导入与展示。
2. 变形的核心算法和参数还需要实验确定。
