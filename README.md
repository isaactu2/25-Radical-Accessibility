# Accessible Design Scripts Repository  
**Grasshopper + Python Toolkit for Architectural Accessibility**

This repository serves as an open-source collection of **parametric design scripts** and **computational workflows** focused on improving accessibility in architectural design studios.  
Each project explores how algorithmic design methods — primarily through Grasshopper and Python — can address spatial accessibility, inclusive circulation, and adaptable environments across architectural typologies.

---

## Overview

The goal of this repository is to make architectural scripting **more accessible**, both conceptually and technically.  
These tools are built to help designers:
- Generate adjustable geometries that meet accessibility standards  
- Automate dimensioning, slope, and path-of-travel analysis  
- Share design logic transparently for iteration and critique  
- Encourage collaborative design learning between studios and disciplines  

While the **Three-Rectangle Gable House Prototypes** are the current focus, the repository is meant to grow into a **hub for accessibility-driven design scripts** for future projects.

---

## Repository Structure

Each project or iteration is organized into its own **branch**.  
Branches act as distinct development phases, so collaborators can easily trace design evolution or fork a specific version.

**Typical folder structure inside each branch:**

```
├── scripts/             # Grasshopper (.ghx) and Python (.py) files
├── exports/             # Renderings, analysis diagrams, and results
├── docs/                # Design notes, iteration logs, and diagrams
└── README.md            # Description of that project or iteration
```

### Branch Naming Conventions
- `iteration-01` → Base logic and geometry setup  
- `iteration-02` → Integrated slope/ramp accessibility  
- `iteration-03` → Structural optimization and parametric facade refinement  

Each branch contains its own documentation, making it easier to compare and annotate changes between versions.

---

## Grasshopper + Python Integration

The **Grasshopper** definitions (`.ghx`) and **Python** scripts (`.py`) in this repository work together to:
- Automate parametric accessibility checks  
- Enable customizable slope and circulation logic  
- Allow scripting-based control for geometry generation  
- Provide clarity and repeatability for design iteration  

The scripts are written for **educational and collaborative purposes**, promoting accessibility both in design outcome and in the process of creating digital models.

---

## Viewing and Contributing

1. **Clone or fork** this repository:
   ```bash
   git clone https://github.com/<username>/accessible-design-scripts.git
   ```
2. **Switch to a project or iteration branch:**
   ```bash
   git checkout iteration-02
   ```
3. Open `.ghx` files in **Grasshopper (Rhino 7 or later)** and `.py` files in **GhPython** or a compatible IDE.

### Contribution Guidelines
- Use Pull Requests to propose edits or new iterations.  
- Annotate parameter logic directly within the `.ghx` file or in markdown.  
- Use **Issues** for bug reports, feature ideas, or discussion on accessibility concepts.  
- Label contributions with tags such as `enhancement`, `iteration`, or `accessibility-study`.

---

## Documentation

The **main branch** contains:
- `docs/iterations.md` — overview of all project phases  
- `docs/references.md` — design standards, accessibility references, and frameworks  
- `README.md` — this file, describing repository purpose and structure  

Each branch adds its own local `README.md` to describe project goals, key parameters, and accessibility outcomes.


---

## Keywords

`#Architecture` `#Accessibility` `#ParametricDesign` `#Grasshopper` `#Rhino` `#Python` `#InclusiveDesign` `#ComputationalDesign`
