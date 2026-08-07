# Academic Work Template Documentation

This repository serves as a template for academic projects that combine LaTeX documentation, Jupyter notebooks for data analysis and visualization, and Python source code for computational implementations.

## Overview

This template provides a structured approach to organizing academic work, ensuring reproducibility, clear documentation, and proper separation of concerns between different types of content.

## Project Structure

### High-Level Organization

```mermaid
graph TD
    A[Academic Work Template] --> B[docs/]
    A --> C[notebooks/]
    A --> D[src/]
    A --> E[scripts/]
    A --> F[Configuration Files]
    
    B --> B1[LaTeX Documents]
    C --> C1[Jupyter Notebooks]
    D --> D1[Source Code]
    E --> E1[Utility Scripts]
    F --> F1[Makefile, requirements.txt, etc.]
    
    style A fill:#e1f5fe
    style B fill:#f3e5f5
    style C fill:#e8f5e8
    style D fill:#fff3e0
    style E fill:#fce4ec
    style F fill:#f1f8e9
```

### Detailed Folder Structure

```mermaid
graph LR
    subgraph "Repository Root"
        A[📁 latex-template/]
    end
    
    subgraph "Documentation"
        B[📁 docs/]
        B --> B1[📁 report/]
        B1 --> B2[📄 report_main.tex]
        B1 --> B3[📄 references.bib]
        B1 --> B4[📄 plain.bst]
        B1 --> B5[📄 Makefile]
        B1 --> B6[📁 build/]
        B6 --> B7[Generated Files]
    end
    
    subgraph "Analysis & Visualization"
        C[📁 notebooks/]
        C --> C1[📓 *.ipynb files]
        C --> C2[Data Analysis]
        C --> C3[Visualizations]
    end
    
    subgraph "Source Code"
        D[📁 src/]
        D --> D1[📁 classes/]
        D --> D2[📁 plots/]
        D1 --> D3[🐍 Python Classes]
        D2 --> D4[📊 Plotting Code]
    end
    
    subgraph "Utilities"
        E[📁 scripts/]
        E --> E1[⚙️ Utility Scripts]
    end
    
    subgraph "Configuration"
        F[📄 Makefile]
        G[📄 requirements.txt]
        H[📄 README.md]
        I[📄 setup-hooks.sh]
    end
    
    A --> B
    A --> C
    A --> D
    A --> E
    A --> F
    A --> G
    A --> H
    A --> I
    
    style A fill:#e3f2fd
    style B fill:#f3e5f5
    style C fill:#e8f5e8
    style D fill:#fff3e0
    style E fill:#fce4ec
```

## Folder Descriptions

### 📁 `docs/` - LaTeX Documents

The `docs` folder contains all LaTeX-related files for generating academic reports, papers, and presentations.

```mermaid
graph TD
    A[docs/] --> B[report/]
    B --> C[📄 report_main.tex]
    B --> D[📄 references.bib]
    B --> E[📄 plain.bst]
    B --> F[📄 Makefile]
    B --> G[📁 build/]
    
    C --> C1[Main LaTeX Document]
    D --> D1[Bibliography Database]
    E --> E1[Bibliography Style]
    F --> F1[Build Automation]
    G --> G1[Compiled Output]
    
    style A fill:#f3e5f5
    style B fill:#f8bbd9
    style C fill:#e1bee7
    style D fill:#e1bee7
    style E fill:#e1bee7
    style F fill:#e1bee7
    style G fill:#e1bee7
```

**Contents:**
- `report_main.tex`: Main LaTeX document file
- `references.bib`: Bibliography database in BibTeX format
- `plain.bst`: Bibliography style file
- `Makefile`: Build automation for LaTeX compilation
- `build/`: Directory containing generated files (PDFs, auxiliary files)

**Purpose:**
- Academic report writing
- Paper drafts and manuscripts
- Thesis chapters
- Presentation slides (if using Beamer)

### 📁 `notebooks/` - Jupyter Notebooks

The `notebooks` folder is designed for interactive data analysis, exploration, and visualization.

```mermaid
graph TD
    A[notebooks/] --> B[Data Analysis]
    A --> C[Exploratory Analysis]
    A --> D[Visualizations]
    A --> E[Prototyping]
    A --> F[Results Documentation]
    
    B --> B1[📓 data_analysis.ipynb]
    C --> C1[📓 exploration.ipynb]
    D --> D1[📓 visualization.ipynb]
    E --> E1[📓 prototype.ipynb]
    F --> F1[📓 results.ipynb]
    
    style A fill:#e8f5e8
    style B fill:#c8e6c9
    style C fill:#c8e6c9
    style D fill:#c8e6c9
    style E fill:#c8e6c9
    style F fill:#c8e6c9
```

**Typical Contents:**
- Data loading and preprocessing notebooks
- Exploratory data analysis (EDA)
- Algorithm implementation and testing
- Results visualization and plotting
- Statistical analysis and hypothesis testing

**Purpose:**
- Interactive development and testing
- Data exploration and visualization
- Rapid prototyping of algorithms
- Documentation of analysis workflow
- Generation of figures for reports

### 📁 `src/` - Source Code

The `src` folder contains the main Python source code organized into logical modules.

```mermaid
graph TD
    A[src/] --> B[classes/]
    A --> C[plots/]
    
    B --> B1[📄 __init__.py]
    B --> B2[🐍 Astar.py]
    B --> B3[🐍 Bfs.py]
    B --> B4[🐍 Dfs.py]
    B --> B5[🐍 Dijkstra.py]
    B --> B6[🐍 GraphPlotter.py]
    B --> B7[🐍 Problem.py]
    B --> B8[🐍 TAD.py]
    
    C --> C1[📓 grafs.ipynb]
    C --> C2[📁 img/]
    
    B1 --> B9[Package Initialization]
    B2 --> B10[A* Algorithm Implementation]
    B3 --> B11[Breadth-First Search]
    B4 --> B12[Depth-First Search]
    B5 --> B13[Dijkstra Algorithm]
    B6 --> B14[Graph Visualization]
    B7 --> B15[Problem Definition]
    B8 --> B16[Abstract Data Types]
    
    C1 --> C3[Graph Plotting Notebook]
    C2 --> C4[Generated Images]
    
    style A fill:#fff3e0
    style B fill:#ffcc80
    style C fill:#ffcc80
```

**Structure:**
- `classes/`: Core Python classes and algorithms
  - Algorithm implementations (A*, BFS, DFS, Dijkstra)
  - Problem definition classes
  - Graph plotting utilities
  - Abstract data types
- `plots/`: Plotting and visualization code
  - Jupyter notebooks for plot generation
  - Generated images and figures

**Purpose:**
- Reusable algorithm implementations
- Core computational logic
- Data structures and classes
- Plotting and visualization utilities
- Modular code organization

### 📁 `scripts/` - Utility Scripts

The `scripts` folder contains utility scripts for project automation and maintenance.

```mermaid
graph TD
    A[scripts/] --> B[🔧 pre-push]
    A --> C[⚙️ Other Utilities]
    
    B --> B1[Git Hook Script]
    C --> C1[Build Scripts]
    C --> C2[Testing Scripts]
    C --> C3[Deployment Scripts]
    
    style A fill:#fce4ec
    style B fill:#f8bbd9
    style C fill:#f8bbd9
```

**Purpose:**
- Git hooks for quality control
- Build and deployment automation
- Testing and validation scripts
- Project maintenance utilities

## Workflow Integration

```mermaid
graph TD
    A[Research Phase] --> B[notebooks/]
    B --> C[Development Phase]
    C --> D[src/]
    D --> E[Documentation Phase]
    E --> F[docs/]
    F --> G[Publication]
    
    B --> B1[Data Analysis]
    B --> B2[Prototyping]
    B --> B3[Visualization]
    
    D --> D1[Implementation]
    D --> D2[Testing]
    D --> D3[Refactoring]
    
    F --> F1[Report Writing]
    F --> F2[Bibliography]
    F --> F3[PDF Generation]
    
    style A fill:#e1f5fe
    style B fill:#e8f5e8
    style C fill:#e1f5fe
    style D fill:#fff3e0
    style E fill:#e1f5fe
    style F fill:#f3e5f5
    style G fill:#e8eaf6
```

## Best Practices

### 1. Development Workflow
1. **Exploration**: Start with notebooks for data exploration and algorithm prototyping
2. **Implementation**: Move stable code to `src/` as reusable modules
3. **Documentation**: Document findings and methodology in LaTeX reports
4. **Integration**: Use Makefile for automated builds and testing

### 2. Code Organization
- Keep notebooks focused on specific analysis tasks
- Extract reusable functions and classes to `src/`
- Maintain clear separation between exploration and production code
- Use proper Python package structure with `__init__.py` files

### 3. Documentation Standards
- Write comprehensive LaTeX documentation in `docs/`
- Include proper citations using BibTeX
- Generate figures in notebooks and reference them in reports
- Maintain version control for all components

### 4. Environment Management
- Use `requirements.txt` for Python dependencies
- Include LaTeX compilation instructions in Makefiles
- Set up Git hooks for code quality checks
- Use virtual environments for Python development

## Getting Started

1. **Clone the template**: Start with this repository structure
2. **Set up environment**: Install Python dependencies and LaTeX distribution
3. **Configure tools**: Set up Git hooks and development environment
4. **Start development**: Begin with notebooks for exploration
5. **Build documentation**: Use Makefiles for automated compilation

This template ensures a professional, organized approach to academic work that promotes reproducibility and clear documentation standards.