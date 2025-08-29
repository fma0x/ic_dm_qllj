# Implementação do potencial q-Lennard-Jones no pacote computacional Lammps  

1. Palavras-chaves: Interações Intermoleculares, q-Lennard-Jones, Dinâmica Molecular,C++.  

## Introdução e Justificativa  

A dinâmica molecular (DM) é uma das técnicas computacionais mais poderosas para
investigar propriedades termodinâmicas, estruturais e dinâmicas de sistemas de partículas
interagentes, com aplicações que vão desde a ciência dos materiais até a biofísica [1].
Tradicionalmente, o potencial de Lennard-Jones (LJ) tem sido amplamente utilizado para
descrever interações van der Waals em sistemas simples, como gases nobres e líquidos
[2]. No entanto, em condições extremas — como altas pressões ou em materiais com
interações não convencionais —, o LJ clássico apresenta limitações que demandam a
adoção de potenciais modificados [3,4]. Nesse contexto, o potencial q-Lennard-Jones (qLJ)
surge como uma generalização promissora, capaz de capturar comportamentos não
triviais em sólidos e fluidos complexos, ajustando-se a cenários onde o LJ padrão falha
[5]. 

A implementação de potenciais personalizados em ferramentas de alta performance,
como o LAMMPS (Large-scale Atomic/Molecular Massively Parallel Simulator), é
essencial para simulações realistas em larga escala [6]. A integração desse potencial em
pacotes otimizados como o LAMMPS permitirá estudos sistemáticos em sistemas com
milhões de partículas, aproveitando paralelismo e otimizações de hardware [7]. Essa
transição de protótipos acadêmicos para ferramentas profissionais é um passo crítico para
validar a utilidade do q-LJ em pesquisas aplicadas, como o desenvolvimento de materiais
sob pressão ou a modelagem de nanoestruturas [8].  

A motivação deste projeto alinha-se a três lacunas identificadas na literatura: (i) a ausência
da implementação do potencial q-Lennard-Jones (q-LJ) em pacotes consolidados de DM
[9]; (ii) a necessidade de validação rigorosa desse potencial em sistemas reais,
comparando-os com dados experimentais [10]; e (iii) a carência de documentação
acessível para pesquisadores que desejem utilizar tais potenciais em seus trabalhos [11].
Além disso, ao abordar essas questões, este projeto contribuirá para o avanço
metodológico da área [12]. 

Por fim, a relevância desta proposta estende-se além do âmbito acadêmico. A capacidade
de simular materiais sob condições extremas com precisão tem implicações diretas em
setores estratégicos, como energia (ex.: armazenamento de hidrogênio) [13] e
nanotecnologia (ex.: design de materiais 2D) [14]. Ao disponibilizar o q-LJ no LAMMPS
— um software amplamente adotado pela comunidade científica [15] —, este projeto visa
criar um legado técnico que facilitará investigações futuras, promovendo colaborações
interdisciplinares, principalmente na investigação de novos materiais [16].  

## Objetivo Geral  

Implementar e validar o potencial q-Lennard-Jones no LAMMPS, permitindo simulações
eficientes de materiais sob condições extremas.
Objetivos Específicos
1. Desenvolver um conjunto de alterações (patch) ou plugin para o LAMMPS que
incorpore o potencial q-LJ.
2. Reproduzir resultados conhecidos do LJ clássico (caso limite).
3. Comparar com dados experimentais (ex.: curvas pressão-volume para metais).
4. Avaliar escalabilidade em sistemas grandes (N > 10⁶ partículas) usando
paralelismo.
5. Documentar e disponibilizar o código no repositório oficial do LAMMPS ou em
plataformas como GitHub.  

## Materiais e Métodos  

Para implementar o potencial q-LJ no LAMMPS, será desenvolvido um pair style
personalizado em C++, seguindo a estrutura modular do código-fonte do LAMMPS [17].
O potencial será integrado como uma nova classe no diretório `src/MOLECULE/`, com
parâmetros ajustáveis via arquivos de entrada. A validação inicial consistirá em reproduzir
resultados do LJ clássico (n = 1, q = 3) para sistemas modelo (e.g., argônio), comparando
propriedades termodinâmicas (energia, pressão) com dados da literatura [18]. Em
seguida, simulações de metais (Fe, Cu) sob alta pressão serão realizadas, utilizando
condições periódicas e termostatos/barostatos acoplados (NPT/NVT), com comparação
direta a curvas pressão-volume experimentais [19].  

A avaliação de desempenho será feita em sistemas com N > 10⁶ partículas, medindo
tempos de execução em diferentes configurações de hardware (CPU/GPU) e comparando
com o LJ nativo do LAMMPS [20]. O código será documentado no formato padrão do
LAMMPS (Doxygen) e disponibilizado em um repositório público, incluindo exemplos
de scripts e tutoriais para usuários finais. Todas as simulações utilizarão o método de
Verlet para integração temporal [21], garantindo estabilidade numérica.  

## Resultados Esperados  

Como resultados concretos deste projeto, espera-se: (1) um pair style funcional para o
potencial q-Lennard-Jones plenamente integrado ao LAMMPS, validado através da
reprodução do limite clássico (LJ padrão) e de simulações de metais (Fe, Cu) sob alta
pressão com erro inferior a 5% em relação a dados experimentais; (2) testes de
desempenho demonstrando eficiência computacional comparável ao LJ nativo em
sistemas com >1 milhão de partículas em arquiteturas paralelas (CPU/GPU); e (3)
documentação completa no formato LAMMPS, incluindo tutoriais e exemplos de 
arquivos de entrada em repositório público, garantindo acessibilidade e reprodutibilidade
pela comunidade científica.  

## Referências  

[1] Allen, M.P.; Tildesley, D.J. Computer Simulation of Liquids. Oxford University Press,
2017.
[2] Rapaport, D.C. The Art of Molecular Dynamics Simulation. Cambridge University
Press, 2004.
[3] Alves, Á.S. et al. A Proposal of an Isothermal q-EoS for Solids at High Pressures.
Sitientibus Ser. Cien. Fis. 20, 1 (2024).
[4] Á. S. Alves, F. A. Farias, R. N. dos Santos, A Three Parameter q-EoS for Solids under
High Pressures. Sitientibus Ser, Cien. Fis. 20, 1 (2024).
[5] Tsallis, C. Introduction to Nonextensive Statistical Mechanics. Springer, 2009.
[6] Plimpton, S. Fast Parallel Algorithms for Short-Range Molecular Dynamics. J. Comp.
Phys. 117, 1-19 (1995).
[7] Thompson, A.P. et al. LAMMPS - a Flexible Simulation Tool for Particle-Based
Materials Modeling. Comp. Phys. Comm. 271, 108171 (2022).
[8] Frenkel, D.; Smit, B. Understanding Molecular Simulation. Academic Press, 2023.
[9] LAMMPS Documentation. Adding a New Pair Style. 2023.
[10] Dewaele, A. et al. High-Pressure-High-Temperature Equations of State for Fe and
MgO. Phys. Rev. B 82, 104106 (2010).
[11] Prlić, A. et al. BioJava: An Open-Source Framework for Bioinformatics.
Bioinformatics 28, 2693-2695 (2012).
[12] Virtanen, P. et al. SciPy 1.0: Fundamental Algorithms for Scientific Computing in
Python. Nature Methods 17, 261-272 (2020).
[13] McMahon, M.I. High-Pressure Physics: The Diamond Age. Nat. Mater. 5, 763-764
(2006).
[14] Geim, A.K. Graphene: Status and Prospects. Science 324, 1530-1534 (2009).
[15] Phillips, J.C. et al. Scalable Molecular Dynamics with NAMD. J. Comp. Chem. 26,
1781-1802 (2005).
[16] Harris, C.R. et al. Array Programming with NumPy. Nature 585, 357-362 (2020).
[17] LAMMPS Developer Guide. Adding a new pair style. Disponível em:
https://docs.lammps.org/Developer_plugins.html. Acessado em 05/04/2025.
[18] Allen, M. P.; Tildesley, D. J. Computer Simulation of Liquids. Oxford, 2017.
[19] Dewaele, A. et al. High-Pressure Equations of State for Fe and Cu. Phys. Rev. B 82,
104106 (2010). 
[20] Thompson, A. P. et al. LAMMPS - a Flexible Simulation Tool. Comp. Phys. Comm.
271, 108171 (2022).
[21] Frenkel, D.; Smit, B. Understanding Molecular Simulation. Academic Press, 2023. 
