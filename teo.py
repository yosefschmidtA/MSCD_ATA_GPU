import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.interpolate import griddata

def calcular_r_factor(df):
    resultados = []
    angulos = []
    resultado_final4 = 0
    n_total = 0

    angulos_unicos = df['Theta'].unique()

    # Acumuladores para o cálculo do Rtotal
    todas_intensidades_exp = []
    todas_intensidades_calc = []

    for theta in angulos_unicos:
        subset = df[df['Theta'] == theta]
        phi = subset['Phi'].values
        intensidade_calculada = np.round(subset['intensitycal'].values, 10)
        intensidade_experimental = np.round(subset['intensityexp'].values, 10)

        # Acumular todas as intensidades para o cálculo do Rtotal
        todas_intensidades_exp.extend(intensidade_experimental)
        todas_intensidades_calc.extend(intensidade_calculada)

        # Normalização por theta
        somatorio1 = np.round(np.sum(np.abs(intensidade_experimental)),10)
        somatorio2 = np.round(np.sum(np.abs(intensidade_calculada)),10)

        if somatorio1 == 0 or somatorio2 == 0:
            continue

        nova_coluna_exp = np.round((intensidade_experimental / somatorio1),10)
        nova_coluna_calc = np.round((intensidade_calculada / somatorio2),10)

        resultado_final1 = np.round(np.sum((nova_coluna_exp - nova_coluna_calc) ** 2),10)
        resultado_final2 = np.round(np.sum(nova_coluna_exp ** 2 + nova_coluna_calc ** 2),10)

        resultado_final3 = np.round((resultado_final1 / resultado_final2),10)
        resultado_final4 += resultado_final3

        resultados.append(resultado_final3)
        angulos.append(theta)

        n_total += 1

    # Cálculo do R-factor médio
    r_factor_medio = np.round((resultado_final4 / n_total),10) if n_total > 0 else None

    # Cálculo do Rtotal usando todas as intensidades acumuladas
    todas_intensidades_exp = np.round(np.array(todas_intensidades_exp),10)
    todas_intensidades_calc = np.round(np.array(todas_intensidades_calc),10)

    somatorio1_total = np.round(np.sum(np.abs(todas_intensidades_exp)),10)
    somatorio2_total = np.round(np.sum(np.abs(todas_intensidades_calc)),10)

    if somatorio1_total > 0 and somatorio2_total > 0:
        nova_coluna_exp_total = np.round((todas_intensidades_exp / somatorio1_total),10)
        nova_coluna_calc_total = np.round((todas_intensidades_calc / somatorio2_total),10)

        resultado_final1_total = np.round(np.sum((nova_coluna_exp_total - nova_coluna_calc_total) ** 2),10)
        resultado_final2_total = np.round(np.sum(nova_coluna_exp_total ** 2 + nova_coluna_calc_total ** 2),10)

        r_factor_total = np.round((resultado_final1_total / resultado_final2_total),10)
    else:
        r_factor_total = None

    return resultados, angulos, r_factor_medio, r_factor_total

def process_file(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()

    data = []
    theta_value = None
    lista_temporaria = []

    for i in range(26, len(lines)):
        line = lines[i].strip()

        if "fitted parameters (" in line:
            break

        if line:
            parts = line.split()

            if len(parts) == 7:
                # Encontrou um novo Theta, processa os dados acumulados
                if lista_temporaria:
                    intensi2_values = [item[1] for item in lista_temporaria]
                    media_intensi2 = sum(intensi2_values) / len(intensi2_values)

                    for phi, intensi1, intensi2, intensityexp in lista_temporaria:
                        intensitycal = (intensi1 - media_intensi2) / media_intensi2
                        data.append([phi, intensitycal, theta_value, intensityexp, True])

                    lista_temporaria = []

                theta_value = float(parts[3])

            elif len(parts) == 5 and theta_value is not None:
                phi = float(parts[0])
                intensi1 = float(parts[1])
                intensi2 = float(parts[2])
                intensityexp = float(parts[4])
                lista_temporaria.append((phi, intensi1, intensi2, intensityexp))

    # Depois do loop, processa o último conjunto também
    if lista_temporaria:
        intensi2_values = [item[1] for item in lista_temporaria]
        media_intensi2 = sum(intensi2_values) / len(intensi2_values)

        for phi, intensi1, intensi2, intensityexp in lista_temporaria:
            intensitycal = (intensi1 - media_intensi2) / media_intensi2
            data.append([phi, intensitycal, theta_value, intensityexp, True])

    df = pd.DataFrame(data, columns=['Phi', 'intensitycal', 'Theta', 'intensityexp', 'IsOriginal'])
    # 1. CÁLCULO DO R-FACTOR (Nos dados originais)
    resultados, angulos, r_factor_medio, r_factor_total = calcular_r_factor(df)
    for theta, r_factor in zip(angulos, resultados):
        print(f'Theta: {theta}, R-factor: {r_factor}')

    print(f'R-factor médio: {r_factor_medio}')
    print(f'R-factor total: {r_factor_total}')

    # -------------------------------------------------------------------------
    # 2. LÓGICA DE SIMETRIA UNIVERSAL (INSERIDA AQUI)
    # -------------------------------------------------------------------------
    # Expande os dados (ex: 0-90) para 360 ANTES de criar os pontos fantasmas.
    phi_interval = df['Phi'].max() - df['Phi'].min()
    step = None
    
    if 170 < phi_interval < 190: 
        step = 180
    elif 80 < phi_interval < 135:
        step = 120 if phi_interval > 105 else 90

    if step is not None:
        dfs_list = [df]
        current_shift = step
        while current_shift < 360:
            df_new = df.copy()
            # Módulo 360 para garantir que 485 vire 125, etc.
            df_new['Phi'] = (df_new['Phi'] + current_shift) % 360
            dfs_list.append(df_new)
            current_shift += step
        # Atualiza o df com a versão expandida
        df = pd.concat(dfs_list, ignore_index=True)
    # -------------------------------------------------------------------------

    # 3. LÓGICA DE PONTOS FANTASMAS (GHOST POINTS) - A SUA LÓGICA
    # Agora funciona perfeitamente pois o 'df' já tem dados perto de 360
    df_plot = df.copy()
    margin = 20 # Margem de 20 graus para cada lado
    
    # Pega o final (ex: 340-360) e joga para trás (-20 a 0)
    subset_end = df_plot[df_plot['Phi'] > (360 - margin)].copy()
    subset_end['Phi'] = subset_end['Phi'] - 360
    
    # Pega o início (ex: 0-20) e joga para frente (360 a 380)
    subset_start = df_plot[df_plot['Phi'] < margin].copy()
    subset_start['Phi'] = subset_start['Phi'] + 360
    
    # Junta tudo (df_plot agora vai de -20 a 380)
    df_plot = pd.concat([df_plot, subset_end, subset_start], ignore_index=True)
    
    # Atualiza o df final para plotagem
    df = df_plot.copy()
    
    return df, r_factor_total  # Retorna r_factor_total

def interpolate_data(df, resolution=1000):
    phi = np.radians(df['Phi'])
    theta = np.radians(df['Theta'])
    intensity = df['intensitycal']

    phi_grid = np.linspace(np.min(phi), np.max(phi), resolution)
    theta_grid = np.linspace(np.min(theta), np.max(theta), resolution)

    phi_grid, theta_grid = np.meshgrid(phi_grid, theta_grid)

    intensity_grid = griddata((phi, theta), intensity, (phi_grid, theta_grid), method='cubic')

    return phi_grid, theta_grid, intensity_grid

def plot_polar_interpolated(df, resolution=500, line_position=0.5, my_variable=None, save_path=None):
    # Interpolar os dados
    plt.ion()
    phi_grid, theta_grid, intensity_grid = interpolate_data(df, resolution)

    # Criando o gráfico polar
    fig, ax = plt.subplots(subplot_kw={'projection': 'polar'}, figsize=(10, 8), dpi=100)

    # Plotando a intensidade interpolada
    c = ax.pcolormesh(phi_grid, theta_grid, intensity_grid, shading='gouraud', cmap='afmhot')

    # Definir o limite máximo do eixo theta com base no maior valor de theta nos dados
    max_theta = df['Theta'].max()  # Maior valor de theta presente nos dados
    ax.set_ylim(0, np.radians(max_theta))  # Limitar o eixo radial até o maior valor de theta

    # Adiciona rótulos para os ângulos theta
    theta_ticks = np.linspace(0, max_theta, num=6)
    ax.set_yticks(np.radians(theta_ticks))
    ax.set_yticklabels([f'{int(tick)}°' for tick in theta_ticks])

    # Adicionando a barra de cores
    cbar = fig.colorbar(c, ax=ax, label='', pad=0.08)
    cbar.set_label('', fontsize=22, fontweight='bold')
    cbar.ax.yaxis.set_tick_params(labelsize=30)
    for label in cbar.ax.get_yticklabels():
        label.set_fontweight('bold')
        
    # Adiciona a variável fora do gráfico
    if my_variable is not None:
        fig.text(0.87, 0.03, f'R-factor: {my_variable}', fontsize=34, color='black', ha='right', va='bottom',
                 fontweight='bold')
                 
    Anysotropy = "Anisotropy"
    fig.text(0.94, 0.9,  Anysotropy, fontsize=34, color='black', ha='right', va='bottom',
             fontweight='bold')
             
    # Definir manualmente os ângulos de Phi
    phi_ticks = np.linspace(0, 2 * np.pi, num=9)[:-1]
    phi_labels = [f'{int(np.degrees(tick))}°' for tick in phi_ticks]
    ax.set_xticks(phi_ticks)
    ax.set_xticklabels(phi_labels, fontsize=26, fontweight='bold')

    # Ajustar individualmente os pads
    pad_values = [1, -1, 3, 0, -7, -6, -1, -6]
    for label, pad in zip(ax.get_xticklabels(), pad_values):
        label.set_y(label.get_position()[1] + pad * 0.01)

    ax.tick_params(pad=8)
    plt.yticks(fontsize=0, fontweight='bold')
    plt.draw()

    if save_path:
        plt.savefig(save_path, dpi=200, bbox_inches='tight')

    plt.pause(600)


# Caminho do arquivo
if __name__ == "__main__":
    if len(sys.argv) > 1:
        file_path = sys.argv[1]
        save_path = 'grafico_polar3.png'
        df, r_factor_total = process_file(file_path)
        
        my_variable = "{:.3f}".format(r_factor_total) if r_factor_total is not None else "N/A"
        
        plot_polar_interpolated(df, resolution=500, line_position=0.5, my_variable=my_variable, save_path=save_path)
    else:
        print("Por favor, forneça o caminho do arquivo como argumento.")
