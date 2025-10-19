import numpy as np
import matplotlib.pyplot as plt

# Parâmetros LJ
epsilon = 1.0
sigma = 1.0

# Array de distâncias
r = np.linspace(0.9, 4.0, 500)
r = r[r > 0]

# Potencial LJ 
V_lj = 4 * epsilon * ((sigma/r)**12 - (sigma/r)**6)

# Valores de corte para comparar
cutoffs = [2.0, 2.5, 3.0, 3.5]

plt.figure(figsize=(12, 8))

# Potencial completo
plt.plot(r, V_lj, 'k-', linewidth=3, label='LJ Completo')

# Potenciais truncados
colors = ['red', 'blue', 'green', 'orange']
for i, rc in enumerate(cutoffs):
    V_trunc = np.where(r <= rc, V_lj, 0)
    plt.plot(r, V_trunc, '--', color=colors[i], linewidth=2, 
             label=f'LJ Truncado (r_c = {rc}σ)')
    plt.axvline(rc, color=colors[i], linestyle=':', alpha=0.5)
    
    # Calcular energia no ponto de corte
    V_rc = 4 * epsilon * ((1/rc)**12 - (1/rc)**6)
    plt.plot(rc, V_rc, 'o', color=colors[i], markersize=8)

plt.axhline(0, color='gray', linestyle='-', alpha=0.3)
plt.axvline(2**(1/6), color='purple', linestyle='--', label='Mínimo (2¹/⁶σ)')

plt.xlabel('Distância (r/σ)', fontsize=12)
plt.ylabel('Energia Potencial (V/ε)', fontsize=12)
plt.title('Efeito do Raio de Corte no Potencial de Lennard-Jones', fontsize=14)
plt.ylim(-1.2, 2)
plt.xlim(0.9, 4)
plt.legend()
plt.grid(True, alpha=0.3)

# Inserir tabela de valores
cutoff_data = []
for rc in cutoffs:
    V_rc = 4 * epsilon * ((1/rc)**12 - (1/rc)**6)
    error_pct = abs(V_rc) * 100  # Erro percentual em relação ao mínimo
    cutoff_data.append(f'rc = {rc}σ: V = {V_rc:.4f}ε ({error_pct:.2f}% de -ε)')

plt.text(2.5, 1.0, '\n'.join(cutoff_data), fontsize=10,
         bbox=dict(boxstyle="round,pad=0.5", facecolor="white", alpha=0.8))

plt.tight_layout()
plt.savefig('lj_cutoff_comparison.png', dpi=300, bbox_inches='tight')
plt.show()

# Análise numérica detalhada
print("Análise do Raio de Corte:")
print("="*50)
for rc in cutoffs:
    V_rc = 4 * epsilon * ((1/rc)**12 - (1/rc)**6)
    pairs_ratio = (rc/2.5)**3  # Razão de pares em relação a 2.5σ
    print(f"r_c = {rc}σ:")
    print(f"  • V(r_c) = {V_rc:.4f}ε ({abs(V_rc)*100:.2f}% de ε)")
    print(f"  • Pares relativos a 2.5σ: {pairs_ratio:.2f}x")
    print(f"  • Ganho precisão vs custo: {'BOM' if rc==2.5 else 'RUIM' if rc>3.0 else 'RAZOÁVEL'}")
    print()
