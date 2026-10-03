// To recompile: g++ -Wall -std=c++17 -O2 -o MolVis.exe main.cpp -lfreeglut -lglu32 -lopengl32 -mwindows

#include <GL/freeglut.h>
#include <iostream>
#include <cmath>
#include <string>
#include <vector>
#include <algorithm>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// --- Molecule Data Structures ---
struct Atom3D {
    float x, y, z;
    int atomicNumber; // 6=C, 1=H, 8=O, 17=Cl, 9=F
};

struct Molecule {
    std::string name;
    std::vector<Atom3D> atoms;
};

struct Category {
    std::string name;
    std::vector<Molecule> molecules;
};

std::vector<Category> categories;
int currentCategoryIndex = 0;
int currentMoleculeIndex = 0;

// --- Global State Variables ---
float camRotX = 20.0f;
float camRotY = 30.0f;
float zoomLevel = -10.0f; 
int lastMouseX = -1;
int lastMouseY = -1;
bool isDragging = false;
bool isDropdownOpen = false;

int windowWidth = 1024;
int windowHeight = 768;

// --- IUPAC Nomenclature Helper ---
std::string getPrefix(int n) {
    const char* p[] = {"Meth", "Eth", "Prop", "But", "Pent", "Hex", "Hept", "Oct", "Non", "Dec",
                       "Undec", "Dodec", "Tridec", "Tetradec", "Pentadec"};
    if (n >= 1 && n <= 15) return p[n-1];
    return "Poly";
}

std::string toLower(std::string str) {
    std::transform(str.begin(), str.end(), str.begin(), ::tolower);
    return str;
}

// --- Dynamic Camera ---
void resetCameraForCurrentMolecule() {
    camRotX = 20.0f;
    camRotY = 30.0f;
    
    if (categories.empty() || categories[currentCategoryIndex].molecules.empty()) return;
    const Molecule& m = categories[currentCategoryIndex].molecules[currentMoleculeIndex];
    
    // Calculate bounding box distance to auto-zoom
    float maxDistSq = 0.0f;
    for(const auto& a : m.atoms) {
        float dSq = a.x*a.x + a.y*a.y + a.z*a.z;
        if (dSq > maxDistSq) maxDistSq = dSq;
    }
    
    float maxRadius = sqrt(maxDistSq);
    zoomLevel = - (maxRadius * 2.2f + 5.0f); // Dynamically push camera back
}

// --- Text Rendering ---
void renderText(float x, float y, void *font, const std::string& string) {
    glRasterPos2f(x, y);
    for (const char* c = string.c_str(); *c != '\0'; c++) {
        glutBitmapCharacter(font, *c);
    }
}

// --- Safe Atom Insertion ---
void addAtomUnique(Molecule& m, float x, float y, float z, int atomicNum) {
    for(const auto& at : m.atoms) {
        float dx = at.x - x, dy = at.y - y, dz = at.z - z;
        if(dx*dx + dy*dy + dz*dz < 0.05f) return; // Prevent overlapping geometry
    }
    m.atoms.push_back({x, y, z, atomicNum});
}

// --- Procedural Molecule Generator ---
// bondType: 1=Single, 2=Double(C1=C2), 3=Triple(C1-C2)
// fgType: 0=None, 1=Alcohol, 2=Acid, 3=Ester, 4=Ether, 5=Chloride, 6=Fluoride
Molecule generateChain(std::string name, int numC, int bondType, int fgType) {
    Molecule m; m.name = name;
    
    for (int i = 0; i < numC; i++) {
        float x = i * 1.3f;
        float y = (i % 2 == 0) ? 0.3f : -0.3f;

        // Flatten backbone for double/triple bonds to represent sp2/sp geometry
        if (bondType == 2 && i < 2) y = 0.0f;
        if (bondType == 3 && i < 3) y = 0.0f;

        addAtomUnique(m, x, y, 0.0f, 6);

        // Terminals and Functional Groups at C0
        if (i == 0) {
            if (fgType == 0) { // Hydrocarbons
                if (bondType == 1) addAtomUnique(m, x-0.9f, y, 0.0f, 1);
                if (bondType == 2) { addAtomUnique(m, x-0.7f, y+0.7f, 0.0f, 1); addAtomUnique(m, x-0.7f, y-0.7f, 0.0f, 1); }
                if (bondType == 3) addAtomUnique(m, x-0.9f, y, 0.0f, 1);
            }
            else if (fgType == 1) { // Alcohol
                addAtomUnique(m, x-1.2f, y, 0.0f, 8); 
                addAtomUnique(m, x-1.7f, y+0.5f, 0.0f, 1); 
            }
            else if (fgType == 2) { // Carboxylic Acid
                addAtomUnique(m, x-0.5f, y+1.1f, 0.0f, 8); // =O
                addAtomUnique(m, x-1.2f, y-0.5f, 0.0f, 8); // -O
                addAtomUnique(m, x-1.7f, y-0.2f, 0.0f, 1); // H
            }
            else if (fgType == 3) { // Ester (Methyl)
                addAtomUnique(m, x-0.5f, y+1.1f, 0.0f, 8); // =O
                addAtomUnique(m, x-1.2f, y-0.5f, 0.0f, 8); // -O
                addAtomUnique(m, x-2.2f, y-0.5f, 0.0f, 6); // C (methyl)
                addAtomUnique(m, x-2.5f, y-1.3f, 0.0f, 1); 
                addAtomUnique(m, x-2.5f, y+0.2f, 0.8f, 1); 
                addAtomUnique(m, x-2.5f, y+0.2f, -0.8f, 1); 
            }
            else if (fgType == 4) { // Ether (Methyl)
                addAtomUnique(m, x-1.2f, y, 0.0f, 8); // -O-
                addAtomUnique(m, x-2.4f, y, 0.0f, 6); // C (methyl)
                addAtomUnique(m, x-2.8f, y+0.8f, 0.0f, 1);
                addAtomUnique(m, x-2.8f, y-0.4f, 0.8f, 1);
                addAtomUnique(m, x-2.8f, y-0.4f, -0.8f, 1);
            }
            else if (fgType == 5) addAtomUnique(m, x-1.4f, y, 0.0f, 17); // Chloride
            else if (fgType == 6) addAtomUnique(m, x-1.2f, y, 0.0f, 9);  // Fluoride
        }
        
        // Right Terminal
        if (i == numC - 1) {
            addAtomUnique(m, x+0.9f, y, 0.0f, 1);
        }

        // Z-Axis Hydrogens (filling remaining valences)
        if (bondType == 3 && i < 2) continue; // Triple bond has no Z-axis H
        if (bondType == 2 && i == 0) continue; // C0 of Alkene already has 2 H's
        if (bondType == 2 && i == 1) {
            addAtomUnique(m, x, y+0.9f, 0.0f, 1); // C1 of Alkene has 1 H
            continue;
        }
        if (i == 0 && (fgType == 2 || fgType == 3)) continue; // Acid/Ester C0 has 3 bonds to O

        // Standard sp3 Z-axis Hydrogens
        addAtomUnique(m, x, y + ((i % 2 == 0) ? 0.6f : -0.6f), 0.8f, 1);
        addAtomUnique(m, x, y + ((i % 2 == 0) ? 0.6f : -0.6f), -0.8f, 1);
    }
    
    // Center the molecule on screen
    float cx = (numC > 0) ? (numC - 1) * 1.3f / 2.0f : 0.0f;
    for(auto& a : m.atoms) a.x -= cx;
    
    return m;
}


// --- Database Initialization ---
void initData() {
    // ==========================================
    // CATEGORY 0: Basic Hydrocarbons
    // ==========================================
    Category basicHC; basicHC.name = "Basic Hydrocarbons";
    
    for (int i = 1; i <= 15; i++) basicHC.molecules.push_back(generateChain(getPrefix(i) + "ane", i, 1, 0));
    for (int i = 2; i <= 15; i++) basicHC.molecules.push_back(generateChain("1-" + getPrefix(i) + "ene", i, 2, 0));
    for (int i = 2; i <= 15; i++) basicHC.molecules.push_back(generateChain("1-" + getPrefix(i) + "yne", i, 3, 0));
    categories.push_back(basicHC);

    // ==========================================
    // CATEGORY 1: Functional Groups
    // ==========================================
    Category functional; functional.name = "Functional Groups";
    
    for (int i = 1; i <= 10; i++) functional.molecules.push_back(generateChain(getPrefix(i) + "anol (Alcohol)", i, 1, 1));
    for (int i = 1; i <= 10; i++) functional.molecules.push_back(generateChain(getPrefix(i) + "anoic Acid", i, 1, 2));
    for (int i = 1; i <= 10; i++) functional.molecules.push_back(generateChain("Methyl " + toLower(getPrefix(i)) + "anoate (Ester)", i, 1, 3));
    for (int i = 1; i <= 10; i++) functional.molecules.push_back(generateChain("Methyl " + toLower(getPrefix(i)) + "yl Ether", i, 1, 4));
    for (int i = 1; i <= 10; i++) functional.molecules.push_back(generateChain("1-Chloro" + toLower(getPrefix(i)) + "ane", i, 1, 5));
    for (int i = 1; i <= 10; i++) functional.molecules.push_back(generateChain("1-Fluoro" + toLower(getPrefix(i)) + "ane", i, 1, 6));
    categories.push_back(functional);

    // ==========================================
    // CATEGORY 2: Carbon Nanostructures
    // ==========================================
    Category cns; cns.name = "Carbon Nanostructures";

    // Graphene
    Molecule graphene; graphene.name = "Graphene Sheet";
    float bond = 1.42f;
    for (int r = -2; r <= 2; ++r) {
        for (int c = -3; c <= 3; ++c) {
            float x0 = 1.5f * bond * c, y0 = sqrt(3.0f) * bond * (r + 0.5f * (std::abs(c) % 2));
            addAtomUnique(graphene, x0, y0, 0.0f, 6);
            addAtomUnique(graphene, x0 + bond, y0, 0.0f, 6);
        }
    }
    cns.molecules.push_back(graphene);

    // Nanotube
    Molecule nanotube; nanotube.name = "Carbon Nanotube";
    float rad = 2.5f;
    for (int r = -3; r <= 3; ++r) {
        for (int c = 0; c < 10; ++c) {
            float x0 = 1.5f * bond * c, y0 = sqrt(3.0f) * bond * (r + 0.5f * (c % 2));
            addAtomUnique(nanotube, rad * cos(x0/rad), y0, rad * sin(x0/rad), 6);
            addAtomUnique(nanotube, rad * cos((x0+bond)/rad), y0, rad * sin((x0+bond)/rad), 6);
        }
    }
    cns.molecules.push_back(nanotube);

    // Diamond Lattice
    Molecule diamond; diamond.name = "Diamond Lattice";
    float a = 3.56f; 
    for(int x = -1; x <= 1; x++) {
        for(int y = -1; y <= 1; y++) {
            for(int z = -1; z <= 1; z++) {
                float ox = x*a, oy = y*a, oz = z*a;
                addAtomUnique(diamond, ox, oy, oz, 6);
                addAtomUnique(diamond, ox, oy+a/2, oz+a/2, 6);
                addAtomUnique(diamond, ox+a/2, oy, oz+a/2, 6);
                addAtomUnique(diamond, ox+a/2, oy+a/2, oz, 6);
                addAtomUnique(diamond, ox+a/4, oy+a/4, oz+a/4, 6);
                addAtomUnique(diamond, ox+a/4, oy+3*a/4, oz+3*a/4, 6);
                addAtomUnique(diamond, ox+3*a/4, oy+a/4, oz+3*a/4, 6);
                addAtomUnique(diamond, ox+3*a/4, oy+3*a/4, oz+a/4, 6);
            }
        }
    }
    cns.molecules.push_back(diamond);

    // Fullerene (C60) via mathematically truncated icosahedron
    Molecule fullerene; fullerene.name = "Fullerene (C60 Buckyball)";
    float phi = 1.6180339887f;
    std::vector<Atom3D> icosa;
    for(int i : {-1, 1}) for(int j : {-1, 1}) {
        icosa.push_back({0.0f, (float)i, j*phi, 6});
        icosa.push_back({(float)i, j*phi, 0.0f, 6});
        icosa.push_back({i*phi, 0.0f, (float)j, 6});
    }
    for(size_t i = 0; i < icosa.size(); i++) {
        for(size_t j = 0; j < icosa.size(); j++) {
            if(i == j) continue;
            float dx = icosa[i].x - icosa[j].x, dy = icosa[i].y - icosa[j].y, dz = icosa[i].z - icosa[j].z;
            if(std::abs(dx*dx + dy*dy + dz*dz - 4.0f) < 0.1f) { // Neighbors are distance 2 apart (sq = 4)
                // Truncate edge at 1/3 distance and scale slightly to match C-C visual bonds
                addAtomUnique(fullerene, (2.0f*icosa[i].x + icosa[j].x)*0.6f, 
                                         (2.0f*icosa[i].y + icosa[j].y)*0.6f, 
                                         (2.0f*icosa[i].z + icosa[j].z)*0.6f, 6);
            }
        }
    }
    cns.molecules.push_back(fullerene);
    categories.push_back(cns);
}

// --- Covalent Radius Helper ---
// Returns approximate covalent radius in angstroms to dynamically calculate bonds
float getCovalentRadius(int atomicNumber) {
    switch(atomicNumber) {
        case 1:  return 0.31f; // Hydrogen
        case 6:  return 0.76f; // Carbon
        case 8:  return 0.66f; // Oxygen
        case 9:  return 0.57f; // Fluorine
        case 17: return 1.02f; // Chlorine
        default: return 0.70f; // Fallback
    }
}

// --- Max Valency Helper ---
int getMaxValency(int atomicNumber) {
    switch(atomicNumber) {
        case 1:  return 1; // Hydrogen
        case 6:  return 4; // Carbon
        case 8:  return 2; // Oxygen
        case 9:  return 1; // Fluorine
        case 17: return 1; // Chlorine
        default: return 4; 
    }
}

// --- 3D Rendering Engine (Computational Graph Solver) ---
void renderMolecule(const Molecule& mol) {
    glEnable(GL_LIGHTING);
    glEnable(GL_LIGHT0);
    glEnable(GL_COLOR_MATERIAL);

    // 1. Draw Spheres
    for (const auto& atom : mol.atoms) {
        glPushMatrix();
        glTranslatef(atom.x, atom.y, atom.z);
        switch(atom.atomicNumber) {
            case 6:  glColor3f(0.22f, 0.22f, 0.22f); glutSolidSphere(0.40, 16, 16); break; 
            case 1:  glColor3f(0.95f, 0.95f, 0.95f); glutSolidSphere(0.25, 16, 16); break; 
            case 8:  glColor3f(0.85f, 0.15f, 0.15f); glutSolidSphere(0.38, 16, 16); break; 
            case 17: glColor3f(0.20f, 0.80f, 0.20f); glutSolidSphere(0.45, 16, 16); break; 
            case 9:  glColor3f(0.20f, 0.80f, 0.80f); glutSolidSphere(0.35, 16, 16); break; 
        }
        glPopMatrix();
    }

    // 2. Build the Bonding Graph
    struct PotentialBond { int i, j; float distSq; };
    std::vector<PotentialBond> pbonds;

    // Gather valid close-proximity pairs
    for (size_t i = 0; i < mol.atoms.size(); ++i) {
        for (size_t j = i + 1; j < mol.atoms.size(); ++j) {
            if (mol.atoms[i].atomicNumber == 1 && mol.atoms[j].atomicNumber == 1) continue;

            float dx = mol.atoms[i].x - mol.atoms[j].x;
            float dy = mol.atoms[i].y - mol.atoms[j].y;
            float dz = mol.atoms[i].z - mol.atoms[j].z;
            float distSq = dx*dx + dy*dy + dz*dz;
            
            float r1 = getCovalentRadius(mol.atoms[i].atomicNumber);
            float r2 = getCovalentRadius(mol.atoms[j].atomicNumber);
            float maxBondDist = (r1 + r2) * 1.35f; 
            
            if (distSq > 0.05f && distSq < maxBondDist * maxBondDist) {
                pbonds.push_back({(int)i, (int)j, distSq});
            }
        }
    }

    // Sort by distance to prevent cross-wiring
    std::sort(pbonds.begin(), pbonds.end(), [](const PotentialBond& a, const PotentialBond& b) {
        return a.distSq < b.distSq;
    });

    struct Bond { int i, j, order; };
    std::vector<Bond> finalBonds;
    std::vector<int> valencyCount(mol.atoms.size(), 0);

    // Form initial single bonds
    for (const auto& pb : pbonds) {
        int v1 = getMaxValency(mol.atoms[pb.i].atomicNumber);
        int v2 = getMaxValency(mol.atoms[pb.j].atomicNumber);
        
        // Strict Valency Enforcement
        if (valencyCount[pb.i] < v1 && valencyCount[pb.j] < v2) {
            finalBonds.push_back({pb.i, pb.j, 1});
            valencyCount[pb.i]++;
            valencyCount[pb.j]++;
        }
    }

    // Upgrade to Double/Triple bonds based on remaining valency capacity
    bool upgraded;
    do {
        upgraded = false;
        for (auto& b : finalBonds) {
            if (b.order < 3) {
                int v1 = getMaxValency(mol.atoms[b.i].atomicNumber);
                int v2 = getMaxValency(mol.atoms[b.j].atomicNumber);
                
                // If both connected atoms still have free slots, upgrade their bond!
                if (valencyCount[b.i] < v1 && valencyCount[b.j] < v2) {
                    b.order++;
                    valencyCount[b.i]++;
                    valencyCount[b.j]++;
                    upgraded = true;
                }
            }
        }
    } while (upgraded);

    // 3. Render Bonds with 3D Perpendicular Offsets
    glDisable(GL_LIGHTING); 
    glColor3f(0.5f, 0.5f, 0.5f); 
    glLineWidth(3.0f);
    glBegin(GL_LINES);

    for (const auto& b : finalBonds) {
        const auto& a1 = mol.atoms[b.i];
        const auto& a2 = mol.atoms[b.j];
        
        if (b.order == 1) {
            glVertex3f(a1.x, a1.y, a1.z);
            glVertex3f(a2.x, a2.y, a2.z);
        } else {
            // Calculate Vector Direction
            float dx = a2.x - a1.x, dy = a2.y - a1.y, dz = a2.z - a1.z;
            float len = sqrt(dx*dx + dy*dy + dz*dz);
            float nx = dx/len, ny = dy/len, nz = dz/len;
            
            // Choose an arbitrary 'up' reference vector to cross-multiply
            float rx = 0, ry = 1, rz = 0;
            if (std::abs(ny) > 0.9f) { rx = 1; ry = 0; } // Fallback if pointing straight up
            
            // Cross product gives us the perpendicular offset vector
            float cx = ry*nz - rz*ny;
            float cy = rz*nx - rx*nz;
            float cz = rx*ny - ry*nx;
            
            float clen = sqrt(cx*cx + cy*cy + cz*cz);
            cx /= clen; cy /= clen; cz /= clen;
            
            float offAmt = 0.12f; // Distance between multiple bond lines
            cx *= offAmt; cy *= offAmt; cz *= offAmt;

            if (b.order == 2) {
                // Two parallel shifted lines
                glVertex3f(a1.x + cx, a1.y + cy, a1.z + cz);
                glVertex3f(a2.x + cx, a2.y + cy, a2.z + cz);
                
                glVertex3f(a1.x - cx, a1.y - cy, a1.z - cz);
                glVertex3f(a2.x - cx, a2.y - cy, a2.z - cz);
            } else if (b.order == 3) {
                // One center line, two shifted lines
                glVertex3f(a1.x, a1.y, a1.z);
                glVertex3f(a2.x, a2.y, a2.z);
                
                float cX3 = cx * 1.5f, cY3 = cy * 1.5f, cZ3 = cz * 1.5f;
                glVertex3f(a1.x + cX3, a1.y + cY3, a1.z + cZ3);
                glVertex3f(a2.x + cX3, a2.y + cY3, a2.z + cZ3);
                
                glVertex3f(a1.x - cX3, a1.y - cY3, a1.z - cZ3);
                glVertex3f(a2.x - cX3, a2.y - cY3, a2.z - cZ3);
            }
        }
    }
    glEnd();
}
// --- Minimalist Overlay ---
void displayUI() {
    glDisable(GL_DEPTH_TEST);
    glMatrixMode(GL_PROJECTION); glPushMatrix(); glLoadIdentity();
    gluOrtho2D(0.0, windowWidth, 0.0, windowHeight);
    glMatrixMode(GL_MODELVIEW); glPushMatrix(); glLoadIdentity();

    // Stats Text
    glColor3f(1.0f, 1.0f, 1.0f);
    renderText(20.0f, windowHeight - 35.0f, GLUT_BITMAP_HELVETICA_18, categories[currentCategoryIndex].molecules[currentMoleculeIndex].name);
    
    glColor3f(0.5f, 0.5f, 0.5f);
    std::string countStr = "Structure " + std::to_string(currentMoleculeIndex + 1) + " of " + std::to_string(categories[currentCategoryIndex].molecules.size());
    renderText(20.0f, windowHeight - 55.0f, GLUT_BITMAP_HELVETICA_12, countStr);
    renderText(20.0f, 20.0f, GLUT_BITMAP_HELVETICA_12, "Use <- / -> arrow keys to toggle structures | Drag mouse to spin | Scroll wheel to zoom");

    // Dropdown UI
    float mx1 = windowWidth - 280.0f, mx2 = windowWidth - 20.0f, my1 = windowHeight - 45.0f, my2 = windowHeight - 15.0f;
    glColor3f(0.12f, 0.12f, 0.15f); glBegin(GL_QUADS); glVertex2f(mx1, my1); glVertex2f(mx2, my1); glVertex2f(mx2, my2); glVertex2f(mx1, my2); glEnd();
    glColor3f(0.35f, 0.35f, 0.4f); glLineWidth(1.0f); glBegin(GL_LINE_LOOP); glVertex2f(mx1, my1); glVertex2f(mx2, my1); glVertex2f(mx2, my2); glVertex2f(mx1, my2); glEnd();
    glColor3f(0.9f, 0.9f, 0.9f); renderText(mx1 + 15.0f, my1 + 8.0f, GLUT_BITMAP_HELVETICA_12, categories[currentCategoryIndex].name + "  [V]");

    if (isDropdownOpen) {
        for (size_t i = 0; i < categories.size(); ++i) {
            float boxY1 = my1 - 30.0f * (i + 1), boxY2 = my1 - 30.0f * i;
            if ((int)i == currentCategoryIndex) glColor3f(0.2f, 0.25f, 0.35f); else glColor3f(0.08f, 0.08f, 0.1f);
            glBegin(GL_QUADS); glVertex2f(mx1, boxY1); glVertex2f(mx2, boxY1); glVertex2f(mx2, boxY2); glVertex2f(mx1, boxY2); glEnd();
            glColor3f(0.2f, 0.2f, 0.25f); glBegin(GL_LINE_LOOP); glVertex2f(mx1, boxY1); glVertex2f(mx2, boxY1); glVertex2f(mx2, boxY2); glVertex2f(mx1, boxY2); glEnd();
            glColor3f(0.8f, 0.8f, 0.8f); renderText(mx1 + 15.0f, boxY1 + 8.0f, GLUT_BITMAP_HELVETICA_12, categories[i].name);
        }
    }

    glPopMatrix(); glMatrixMode(GL_PROJECTION); glPopMatrix(); glMatrixMode(GL_MODELVIEW); glEnable(GL_DEPTH_TEST);
}

// --- Main Drawing ---
void display() {
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT); glLoadIdentity();
    glTranslatef(0.0f, 0.0f, zoomLevel); glRotatef(camRotX, 1.0f, 0.0f, 0.0f); glRotatef(camRotY, 0.0f, 1.0f, 0.0f);
    
    if (!categories.empty()) renderMolecule(categories[currentCategoryIndex].molecules[currentMoleculeIndex]);
    
    displayUI(); glutSwapBuffers();
}

// --- Inputs ---
void mouse(int button, int state, int x, int y) {
    float mx = (float)x, my = (float)(windowHeight - y);
    float mx1 = windowWidth - 280.0f, mx2 = windowWidth - 20.0f, my1 = windowHeight - 45.0f, my2 = windowHeight - 15.0f;

    if (button == GLUT_LEFT_BUTTON) {
        if (state == GLUT_DOWN) {
            if (mx >= mx1 && mx <= mx2 && my >= my1 && my <= my2) { isDropdownOpen = !isDropdownOpen; glutPostRedisplay(); return; }
            if (isDropdownOpen && mx >= mx1 && mx <= mx2) {
                for (size_t i = 0; i < categories.size(); ++i) {
                    if (my >= my1 - 30.0f * (i + 1) && my < my1 - 30.0f * i) {
                        currentCategoryIndex = (int)i; currentMoleculeIndex = 0; isDropdownOpen = false;
                        resetCameraForCurrentMolecule(); glutPostRedisplay(); return;
                    }
                }
            }
            isDropdownOpen = false; isDragging = true; lastMouseX = x; lastMouseY = y;
        } else isDragging = false;
    } else if (button == 3 || button == 4) {
        if (state == GLUT_DOWN) { zoomLevel += (button == 3) ? 1.0f : -1.0f; glutPostRedisplay(); }
    }
}

void motion(int x, int y) {
    if (isDragging) {
        camRotY += (x - lastMouseX) * 0.4f; camRotX += (y - lastMouseY) * 0.4f;
        lastMouseX = x; lastMouseY = y; glutPostRedisplay();
    }
}

void special(int key, int x, int y) {
    auto& group = categories[currentCategoryIndex].molecules;
    if (key == GLUT_KEY_RIGHT) currentMoleculeIndex = (currentMoleculeIndex + 1) % group.size();
    else if (key == GLUT_KEY_LEFT) currentMoleculeIndex = (currentMoleculeIndex - 1 + group.size()) % group.size();
    resetCameraForCurrentMolecule(); glutPostRedisplay();
}

void reshape(int w, int h) {
    if (h == 0) h = 1; windowWidth = w; windowHeight = h;
    glViewport(0, 0, w, h); glMatrixMode(GL_PROJECTION); glLoadIdentity();
    gluPerspective(45.0, (float)w / (float)h, 0.1, 100.0); glMatrixMode(GL_MODELVIEW);
}

int main(int argc, char** argv) {
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH | GLUT_MULTISAMPLE); // Multi-sample strictly retained
    glEnable(GLUT_MULTISAMPLE);
    glutInitWindowSize(windowWidth, windowHeight);
    glutCreateWindow("Full Organic Visualizer");
    glEnable(GL_DEPTH_TEST);
    glClearColor(0.04f, 0.04f, 0.05f, 1.0f); 

    GLfloat lightPos[] = { 12.0f, 15.0f, 10.0f, 1.0f }; GLfloat difLight[] = { 0.85f, 0.85f, 0.85f, 1.0f };
    glLightfv(GL_LIGHT0, GL_POSITION, lightPos); glLightfv(GL_LIGHT0, GL_DIFFUSE, difLight);

    initData();
    resetCameraForCurrentMolecule();
    glutDisplayFunc(display);
    glutReshapeFunc(reshape);
    glutMouseFunc(mouse);
    glutMotionFunc(motion);
    glutSpecialFunc(special);
    glutMainLoop();
    return 0;
}