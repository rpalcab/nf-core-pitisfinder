# Usa una imagen base de Miniconda
FROM continuumio/miniconda3

# Establece variables de entorno para Matplotlib y Fontconfig
ENV MPLCONFIGDIR=/tmp/matplotlib_cache
ENV HOME=/home/appuser

# Crea un usuario no root y los directorios necesarios
RUN useradd -m -s /bin/bash appuser && \
    mkdir -p $MPLCONFIGDIR && \
    chown -R appuser:appuser $MPLCONFIGDIR

# Copia el archivo environment.yml al contenedor
COPY environment.yml /tmp/environment.yml

# Instala el entorno desde environment.yml
RUN conda env create -n visualizer -f /tmp/environment.yml

# Activa el entorno de Conda y establece el entorno como predeterminado
RUN echo "source activate visualizer" > ~/.bashrc
ENV PATH=/opt/conda/envs/visualizer/bin:$PATH

# Cambia al usuario no root
USER appuser

# Establece el directorio de trabajo
WORKDIR /home/appuser

# Configura el comando de inicio (opcional, puede ser bash u otro script)
CMD ["bash"]
