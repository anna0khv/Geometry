
///////////////////////////////////////////////////// 9.1

bool pointInRectangle(Point& p, Rectangle& R)
{
	return ((R.sw.x <= p.x) && (p.x <= R.ne.x) && (R.sw.y <= p.y) && (p.y <= R.ne.y));
}

//////////////////////////////////////////////////////// 9.2.1


class Grid {
private:
	int m; // Число ячеек на сторону домена
	double cellSize; // размер стороны ячейки
	List<Point*>*** g; // сетка m x m элементов
	void _Grid(double domainSize, Point s[], int n);
public:
	Grid(double domainSize, Point s[], int n, int m = 10);
	Grid(double domainSize, Point s[], int n, double M = 1.0);
	~Grid(void);
	List<Point*>* rangeQuery(Rectangle& range);
	frind class Quadtree;
};

//////////////////////////////////////////////////////

Grid::Grid(double domainSize, Point s[], int n, int _m) : 9.2.2
	m(_m)

{
	_Grid(domainSize, s, n);
}

//////////////////////////////////////////////////////

void Grid::_Grid(double domainSize, Point s[], int n)
{
	cellSize = domainSize / m;
	g = new (List<Point*>**)[m]; 
	for (int i = 0; i < m; i++)
	{
		g[i] = new (List<Point*>*)[m];
		for (int j = 0; j < m; j++)
			g[i][j] = new List<Point*>;
	}
	for (int i = 0; i < n; i++)
	{
		int a = int(s[i].x / cellSize); 
		int b = int(s[i].y / cellSize);
		g[a][b]->append(new Point(s[i]));
	}
}

///////////////////////////////////////////////////////

Grid::Grid(double domainSize, Point s[], int n, double M) :
	m(int(ceil(sqrt(n / M))))
{
	_Grid(domainSize, s, n);
}

////////////////////////////////////////////////////////

Grid::~Grid(void)
{
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < m; j++) {
			g[i][j]->last();
			while (g[i][j]->lenght() > 0)
				delete g[i][j]->remove();
			delete g[i][j];
		}
		delete g[i];
	}
	delete g;
}

///////////////////////////////////////////////////////////////// 9.2.3

List<Point*> *Grid::rangeQuery(Rectangle &R)
{
	List<Point*>* result = new List<Point*>;
    int ilimit = int(R.ne.x / cellSize);
	int jlimit = int(R.ne.y / cellSize);
	for(int i = int(R.sw.x / cellSize); i <= ilimit; i++) 
		for (int j = int(R.sw.y / cellSize); j <= jlimit; j++) {
			List<Point*>* pts = g[i][j];
			for (pts->first(); !pts->isHead(); pts->next()) {
				Point* p = pts->val();
				if (pointInRectangle(*p, R))
					result->append(p);
			}
		}
	return result;
}

/////////////////////////////////////////////////////////////// 9.3.1

class Quadtree {
private:
	QuadtreeNode* root;
	Rectangle domain;
	QuadtreeNode* buildQuadtree(Grid& G, int M, int D, int level, int, int, int, int);
public:
	Quadtree(Grid& G, int M, int D);
	~Quadtree();
	List<Point*>* rangeQuery(Rectangle &range);
};

///////////////////////////////////////////////////////////////

class QuadtreeNode {
private:
	QuadtreeNode* child[4];
	List<Point*>* pts; // точки в S, если узел вкешний 
	                     // NULL, если узел внутренний
	int size;          // число точек в S, накрываемых узлом
	List<Point*> #rangeQuery(Rectangle &range, Rectangle &quadrant);
	Rectangle quadrant(Rectangle&, int);
	int isExternal();
public:
	QuadtreeNode(List<Point*>*);
	QuadtreeNode(void);
	~QuadtreeNode(void);
	friend class Quadtree;
};

//////////////////////////////////////////////////////////////// 9.3.2

Quadtree::Quadtree(Grid& G, int M, int D)
{
	root = buildQuadtree(G, M, D, 0, 0, G.m - 1, 0, G.m - 1);
	domain = Rectangle(Point(0, 0),
		Point(G.m * G.cellSize, G.m * G.cellSize));
}

///////////////////////////////////////////////////////////////

QuadtreeNode* Quadtree: : buildQuadtree(Grid 6G, int M, int D, int level, int imin, int imax, int jmin, int jmax)
if (imin == imax) (QuadtreeNode * q = new QuadtreeNode(G.g[imin][jmin]);
G.g(imin][jmin] = new List<Point*>;
return g;
) else (
	QuadtreeNode * p = new QuadtreeNode;
int imid = (imin + imax) / 2;
int jmid = (jmin + jmax) / 2;
p->child[0] =
// северо-восток buildQuadtree (G, M, D, level+1, imid+1, imax, jmid+1, jmax) p->child[1] =
// юго-восток
buildguadtree(G, M, D, level + 1, imid + 1, imax, jmin, jmid) p->child[2] =
// юго-запад
;
;
buildguadtree(G, M, D, level + 1, imin, imid, jmin, jmid);
p->child[3] =
// сезеро-запад
buildQuadtree(G, M, D, level + 1, imin, imid, jmid + 1, jmax);
for (int i = 0; i < 4; i++)
	p->size += p->child[i]->size;
if ((p->size <= M) || (level > = D))( //слияние потонков p->pts = new List<Point*>;
	for (i = 0; i < 4; i++) (
		p->pts->append(p->child[i]->pts); delete p->child[i]; p->child[i] = NULL;
)
// завершение слияния потомков
return p:

///////////////////////////////////////////////////

Quadtree::~Quadtree()
{
	delete root;
}

///////////////////////////////////////////////////

QuadtreeNode::QuadtreeNode(List<Point*>* _pts) :
	pts(_pts), size(_pts->length())
{
	for (int i = 0; i < 4; i++)
		child[i] = NULL;
}

///////////////////////////////////////////////////

QuadtreeNode::QuadtreeNode(void) :
	pts(NULL), size(0)
{
	for (int i = 0; i < 4; i++)
		child[i] = NULL;
}

///////////////////////////////////////////////////

QuadtreeNode::~QuadtreeNode()
{
	if (isExternal()) {    		// узел внешний
		pts->last();
		while (pts->length() > 0)
			delete pts->remove();
		delete pts;
	}
	else				// узел внутренний
		for (int i = 0; i < 4; i++)
			delete child[i];
}

/////////////////////////////////////////////////// 9.3.3

List<Point*>* Quadtree::rangeQuery(Ractangle& R)
{
	return root->rangeQuery(R, domain);
}

/////////////////////////////////////////////////// 

List<Point*>* QuadtreeNode::rangeQuery(Rectangle& R, Rectangle& span)
{
	List<Point*>* result = new List<Point*>;
	if (!intersect(R, span))
		return result;
	else if (isExternal()) // узел внешний
		for (pts->first(); !pts->isHead(); pts->next()) {
			Point* p = pts->val();
			if (pointInRectangle(*p, R))
				result->append(p);
		}
		else
			for (int i = 0; i < 4; i++) (
				List<Point*> *1 =
					child[i]->rangeQuery(R, quadrant(span, i));
			result->append(1);
		}
	return result;
}

/////////////////////////////////////////////////// 9.3.4

Rectangle QuadtreeNode::quadrant(Rectangle& s, int i)
{
	Point c = 0.5 * (s.sw + s.ne);
	switch (i) {
	case 0:
		return Rectangle(c, s.ne);
	case 1:
		return Rectangle(Point(c.x, s.sw.y), Point(s.ne.x, c.y));
	case 2:
		return Rectangle(s.sw, c);
	case 3:
		return Rectangle(Point(s.sw.x, c.y), Point(c.x, s.ns.y));
	}
}

//////////////////////////////////////////////////

bool intersect(Rectangle& a, Rectangle& b)
{
	return (overlappingExtent(a, b, X) &&
		overlappingExtent(a, b, Y));
}

//////////////////////////////////////////////////

bool overlappingExtent(Rectangle &a, Rectangle &b, int i)
{
	return ((a.sw[i] <= b.sw[i]) && (b, sw[i] <= a.ne[i])) ||
		((b.sw[i] <= a.sw[i]) && (a.sw[i] <= b.ne[i]));
}

///////////////////////////////////////////////////

int QuadtreeNode::isExternal()
{
	return (pts != NULL);
}

/////////////////////////////////////////////////// 9.4.1

class TwoDTree {
private:
	TwoDTreeNode* root;
	TwoDTreeNode* buildTwoDTree(Point* x[], Point* y[], int n, int cutType);
public:
	TwoDTree(Point p[], int n);
	~TwoDTree(void);
	List<Point*>* rangeQuery(Rectangle& range);
}


/////////////////////////////////////////////

class TwoDTreeNode {
private:
	Point* pnt; // точка, ассоциированная с узлом 
	TwoDTreeNode* lchild; // слева или под линией отсечения 
	TwoDTreeNode* rchild; // справа или над линией стсечения 
	List<Point*> #rangeQuery(Rectangle& range, int cutType);
public:
	TwoDTreeNode(Point*);
	~TwoDTreeNode(void);
	friend class TwoDTree;
};

/////////////////////////////////////////////// 9.4.2

TwoDTree::TwoDTree(Point p[], int n)
{
	Paint** x = new (Point*)[n];
	Point** y = new (Point*)[n];
	for (int i = 0; i < n; i++)
		x[i] = y[i] = new Point(p[i]);
	mergeSort(x, n, leftToRightCmp);
	mergeSort(y, n, bottomToTopCmp);
	root = buildTwoDTree(x, y, n, VERTICAL);
}

//////////////////////////////////////////////

enum (VERTICAL = 0, HORIZONTAL = 1);

/////////////////////////////////////////////

TwoDTreeNode* TwoDTree::buildTwoDTree(Point* x[], Point* y[], int n, int cutType)
{
	if (n == 0)
		return NULL;
	else if (n == 1)
		return new TwoDTreeNode(x[0]);
	int m = n / 2;
	int (*cmp) (Point*, Point*);
	if (cutType == VERTICAL) cmp = leftToRightCmp;
	else cmp = bottomToTopCmp;
	TwoDTreeNode* p = new TwoDTreeNode(x[m]);
	Point** yL = new (Point*)[m];
	Point** yR = new (Point*)[n - m];
	splitPointSet(y, n, x[m], yL, yR, cmp);
	p->lchild = buildTwoDTree(yL, x, m, 1 - cutType);
	p->rchild = buildTwoDTree(yR, x + m + l, n - m - 1, 1 - cutType);
	delete yL;
	delete yR;
	return p;
}

///////////////////////////////////////////////

TwoDTree::~TwoDTree()
{
	delete root;
}

/////////////////////////////////////////////

TwoDTreeNode::TwoDTreeNode(Point* pnt) :
	pnt(_pnt), lchild(NULL), rchild(NULL)
	(
		)

	////////////////////////////////////////////

	TwoDTreeNoda::~TwoDTreeNode()
{
	if (lchild) delete lchild;
	if (rchild) delete rchild;
	delete pnt;
}

/////////////////////////////////////////// 9.4.3

List<Point*>* TwoDTree::rangeQuery(Rectangle& R)
{
	return root->rangeQuery(R, VERTICAL);
}

//////////////////////////////////////////

List<Point*>* TwoDTreeNode::rangeQuery(Rectangle& R, int cutType)
{
	List<Point*>* result = new List<Point*>;
	if (pointInRectangle(*pnt, R))
		result->append(pnt);
	int (*cmp) (Point*, Point*);
	cmp = (cutType == VERTICAL) ? leftToRightCmp : bottomToTopCmp;
	if (lchild && ((*cmp) (&R.sw, pnt) < 0))
		result->append(lchild->rangeQuery(R, 1 - cutType));
	if (rchild && ((*cmp) (&R.ne, pnt) > 0))
		result->append(rchild->rangeQuery(R, 1 - cutType)) :
		return result;
}


///////////////////////////////////////////// 9.4.4.

void splitPointSet(Point* y[], int n, Point* p,
	Point* yL[], Point* yR[],
	int (*cmp)(Point*, Point*))
{
	int lindx = 0, rindx = 0;

	for (int i = 0; i < n; i++) {
		if ((*cmp) (y[i], p) < 0)
			YL[lindx++] = y[i];
		else if ((*cmp) (y[i], p) > 0)
			yR[rindx++] = y[i];
	}
}

/////////////////////////////////////////////////// 9.5.1

class BspTree {
private:
	BspTreeNode* root;
	BspTreeNods* buildBspTree(List<Triangle3D*>*);
public:
	BspTree(Triangle3D* t[], int n);
	~BspTree(void);
	List<Triangle3D*>* visibilitySort(Point3D p);
};

///////////////////////////////////////////////////

class BspTreeNode {
private:
	BspTreeNode* poschild;
	BspTreeNode* negchild;
	Triangle3D* tri;
	BspTreeNode(Triangle3D*);
	~BspTreeNode(void);
	List<Triangle3D*>* visibilitySort(Point3D);
	friend class BspTree;
};

/////////////////////////////////////////////////// 9.5.2

BspTree::BspTree(Triangle3D* t[], int n)
{
	List<Triangle3D*>* this = new List<Triangle3D*>;
	for (int i = 0; i < n; i++)
		this->append(new Triangle3D(*t[i]));
	root = buildBspTree(this);
	delete this;
};

///////////////////////////////////////////////////

BspTreeNode* BspTree::buildBspTree(List<Triangle3D*>* s)
{
	if (s->length() == 0)
		return NULL;
	if (s->length() == 1)
		return new BspTreeNode(s->first());
	List<Triangle3D*>* sP = new List<Triangle3D*>;
	List<Triangle3D*>* sN = new List<Triangle3D*>;
	Triangle3D* p = s->first();
	for (s->next(); !s->isHead(); s->next()) {
		Triangle3D* q = s->val();
		int cl[3];
		for (int i = 0; i < 3; i++)
			cl[i] = (*q)[i].classify(*p);
		if ((cl[0] != NEGATIVE) &&
			(cl[1] != NEGATIVE) &&
			(cl[2] != NEGATIVE))
			sP->append(q);
		else if ((cl[0] != POSITIVE) &&
			(c1[1] != POSITIVE) &&
			(cl[2] != POSITIVE))
			sN->append(q);
		else
			refineList(s, P);
		BspTreeNode* n = new BspTreeNode(p);
		n->poschild = buildBspTree(sP);
		n->negchild = buildBspTree(sN);
		delete sP;
		delete sN;
		return n;
	}
};

///////////////////////////////////////////////////

BspTree::~BspTree()
{
	delete root;
}

///////////////////////////////////////////////////

BspTreeNode::BspTreeNode(Triangle3D* _tri) :
	tri(_tri), poschild(NULL), negchild(NULL)
{
}

///////////////////////////////////////////////////

BspTreeNode::~BspTreeNode()
(
	if (poschild) delete poschild;
if (negchild) delete negchild;
delete tri;
}

/////////////////////////////////////////////////// 9.5.3

List<Triangle3D*>* BspTree::visibilitySort(Point3D p)
{
	return root->visibility Sort(p);
}

///////////////////////////////////////////////////

List<Triangle3D*>* BspTreeNode::visibilitySort(Point3D p)
{
	List<Triangle3D*>* s = new List<Triangle3D*>;
	if (p.classify(*tri) == POSITIVE) {
		if (negchild) s->append(negchild->visibilitySort(p));
		s->append(tri);
		if (poschild) s->append(poschild->visibilitySort(p));
	}
	else {
		if (poschild) s->append(poschild->visibilitySort(p));
		s->append(tri);
		if (neg child) s->append(negchild->visibilitySort(p));
	}
	return s;
}

