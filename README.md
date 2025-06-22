# 🧠 NeuroGraph-Dyscalculia

**Análisis de redes génicas asociadas al fenotipo Discalculia mediante teoría de grafos y enriquecimiento funcional.**

Este proyecto tiene como objetivo identificar y analizar genes relacionados con la discalculia, una alteración del aprendizaje que afecta la capacidad para realizar operaciones matemáticas. Utilizando herramientas de teoría de grafos y análisis de enriquecimiento funcional, se ha construido una red génica que permite explorar relaciones funcionales, estructurales y patológicas entre genes, y su posible impacto en el sistema nervioso.

---

## 📌 Objetivos

- Identificar genes asociados al fenotipo *Discalculia* a través de la base de datos HPO.
- Construir una red de interacciones proteína-proteína con STRING.
- Aplicar técnicas de teoría de grafos para detectar comunidades funcionales.
- Realizar análisis de enriquecimiento funcional con Gene Ontology (GO).

---

## 🧰 Herramientas y Tecnologías

- `R` y librerías: `igraph`, `Linkcomm`, `clusterProfiler`, `org.Hs.eg.db`, `biomaRt`
- Ontologías y recursos:
  - [HPO - Human Phenotype Ontology](https://hpo.jax.org/)
  - [STRING Database](https://string-db.org/)
  - [Gene Ontology (GO)](http://geneontology.org/)

---

## 📊 Metodología

1. **Extracción de genes asociados a Discalculia** desde HPO.
2. **Construcción de red PPI** con genes y sus vecinos usando STRING.
3. **Clustering** mediante LinkComm para detectar comunidades de genes.
4. **Enriquecimiento funcional** con clusterProfiler para identificar funciones biológicas significativas.
5. **Interpretación biológica**: análisis de las conexiones genéticas con procesos neurológicos clave.

---

## Participantes: 
* David Ramírez Arco (Archerd6)
* Alejandro Pascual Mellado (alexpascualm)
