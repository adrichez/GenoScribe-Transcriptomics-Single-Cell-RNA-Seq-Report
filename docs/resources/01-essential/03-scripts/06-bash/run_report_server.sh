#!/bin/bash

set -euo pipefail


#############################################################################################################################################################
# SCRIPT PARA LEVANTAR UN SERVIDOR LOCAL PARA VISUALIZAR EL INFORME HTML CON TODAS LAS FUNCIONALIDADES
#############################################################################################################################################################

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# AJUSTES DE PERSONALIZACIÓN DE LA SALIDA
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#=====================================================================================================================================
# Colores y estilos de texto
#=====================================================================================================================================

# Estilos de texto
BOLD="\033[1m"
UNDERLINE="\033[4m"
ITALIC="\033[3m"


# Color de reinicio
RESET="\033[0m"

# Colores para mensajes
GREY="\033[38;5;245m"
CYAN="\033[0;36m"
YELLOW="\033[38;5;226m"

GREEN="\033[0;32m"
LIME="\033[38;5;118m"
ORANGE="\033[38;5;208m"
RED="\033[0;31m"


# Colores para títulos
BLUE="\033[0;34m"
TEAL="\033[38;5;44m"
PURPLE="\033[38;5;57m"
MAGENTA="\033[38;5;129m"
PINK="\033[38;5;218m"


# Funciones
note() { echo -e "${GREY}${ITALIC}$*${RESET}"; }
info() { echo -e "${CYAN}$*${RESET}"; }
important() { echo -e "${YELLOW}$*${RESET}"; }

success() { echo -e "${GREEN}$*${RESET}"; }
success_b() { echo -e "${LIME}$*${RESET}"; }
warn() { echo -e "${ORANGE}$*${RESET}"; }
error() { echo -e "${RED}$*${RESET}"; }

title1() { echo -e "${BLUE}${BOLD}$*${RESET}"; }
title2() { echo -e "${TEAL}${BOLD}$*${RESET}"; }
title3() { echo -e "${PURPLE}${BOLD}$*${RESET}"; }
title4() { echo -e "${MAGENTA}${BOLD}$*${RESET}"; }
title5() { echo -e "${PINK}${BOLD}$*${RESET}"; }




#=====================================================================================================================================
# Espaciado
#=====================================================================================================================================

TAB_4=$'\033[4G'  # Tabulación a la columna 4
TAB_5=$'\033[5G'  # Tabulación a la columna 5






#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# CODIGO PRINCIPAL
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Configuración de variables robusta
# Cambiar al directorio donde está el script y definir ruta
PATH_SCRIPT_LAUNCH="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$PATH_SCRIPT_LAUNCH"

# Directorio donde está este script
DIR_SCRIPT_LAUNCH="$(cd "$(dirname "$0")" && pwd)"

# Definir la ruta raíz y el nombre del repositorio
PATH_REPORT="$(cd "$PATH_SCRIPT_LAUNCH/../../../.." && pwd)"  # 4 niveles arriba
REPORT_BASENAME="$(basename "$PATH_REPORT")"

# Definir resto de variables
PORT=8000

# URL para abrir en navegador
URL="http://localhost:$PORT/index.html?nocache=$(date +%s)"

echo ""
title1 "###################################################################################################################################"
title1 "🚀${TAB_4}Ejecutanto script para abrir informe HTML a través de un servidor local"
title1 "###################################################################################################################################"
echo ""

title2 "🛠️${TAB_4}Parámetros de configuración:"
title2 "================================"
info "🔹${TAB_4}Directorio del script: $PATH_SCRIPT_LAUNCH"
info "🔹${TAB_4}Ruta raíz del repositorio: $PATH_REPORT"
info "🔹${TAB_4}Nombre del repositorio: $REPORT_BASENAME"
info "🔹${TAB_4}Puerto del servidor: $PORT"

echo ""
title2 "⚙️${TAB_4}Realizando comprobaciones previas..."
title2 "========================================="
# Comprobar Python3
echo "⏳${TAB_4}Comprobando si Python3 está instalado..."
if ! command -v python3 &> /dev/null
then
    error "❌${TAB_4}Python3 no está instalado"
    exit 1
fi
success "✅${TAB_4}Python3 detectado: $(python3 --version)"

echo ""
title2 "⏳${TAB_4}Abriendo navegador web según el sistema operativo..."
title2 "========================================================="
# Abrir navegador según SO
OS="$(uname)"
important "🌐${TAB_4}Intentando abrir el navegador en ${URL}..."
if [[ "$OS" == "Darwin" ]]; then
    echo "💻${TAB_4}Sistema operativo: macOS"
    echo "🚀${TAB_4}Abriendo navegador con..."
    open "$URL"
    success "✅${TAB_4}Navegador abierto correctamente"
elif [[ "$OS" == "Linux" ]]; then
    echo "💻${TAB_4}Sistema operativo: Linux"
    echo "🚀${TAB_4}Abriendo navegador..."
    xdg-open "$URL" 2>/dev/null || echo -e "Abre $URL manualmente"
    success "✅${TAB_4}Navegador abierto correctamente"
else
    warn "💻${TAB_4}Sistema operativo: Desconocido"
    echo -e "Abre $URL manualmente en su navegador"
fi

echo ""
title2 "🚀${TAB_4}Lanzando servidor HTTP con Python3..."
title2 "=========================================="

# Cerrar cualquier proceso Python que esté usando el puerto actual
if lsof -i :$PORT >/dev/null 2>&1; then
    warn "⚠️${TAB_4}Ya hay un servidor ejecutándose en el puerto $PORT. Cerrándolo..."
    kill -9 $(lsof -t -i :$PORT) 2>/dev/null || true
    success "✅${TAB_4}Servidor anterior detenido correctamente."
fi

# Levantar servidor HTTP desde el directorio del informe
important "⏳${TAB_4}Levantando servidor en puerto ${PORT} desde ${PATH_REPORT}..."
important "⏳${TAB_4}Abriendo el navegador automaticamente..."
important "🌐${TAB_4}En el caso de que no se abra automáticamente, acceda a ${URL} manualmente en su navegador"
important "⌨️${TAB_4}Es recomendable pulsar Ctrl+Shift+R o Cmd+Shift+R para forzar la recarga completa de la página y evitar cachés antiguos"
note "💡${TAB_4}Es posible que el navegador no muestre los cambios más recientes debido al caché y en su lugar se muestre un informe antiguo si este se abrió anteriormente mediante la misma URL"
note "🙏${TAB_4}Por este motivo, siempre es recomendable realizar este paso, para estar seguros de que la información está actualizada"
info "⌨️${TAB_4}Presione Ctrl+C para detener el servidor"
echo ""

cd "$PATH_REPORT" || exit

# Ejecutar Python HTTP server
python3 -m http.server "$PORT"

echo ""
