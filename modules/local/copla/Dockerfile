# Usar Miniconda como base
FROM continuumio/miniconda3

# Instalar programas secundarios que dan fallo al construir la imagen

RUN apt-get update && \ 
apt-get install -yq --no-install-recommends \ 
libasound2 libatk1.0-0 libc6 libcairo2 libcups2 libdbus-1-3 \ 
libexpat1 libfontconfig1 libgcc1 libgconf-2-4 libgdk-pixbuf2.0-0 libglib2.0-0 libgtk-3-0 libnspr4 \ 
libpango-1.0-0 libpangocairo-1.0-0 libstdc++6 libx11-6 libx11-xcb1 libxcb1 \ 
libxcursor1 libxdamage1 libxext6 libxfixes3 libxi6 libxrandr2 libxrender1 libxss1 libxtst6 \ 
libnss3

# Definir el directorio de trabajo
WORKDIR /opt/copla

# Descargar COPLA
RUN git clone https://github.com/santirdnd/COPLA .
# Enlace simbólico para que copla no de fallo
RUN ln -s /opt/conda/etc /opt/copla/etc

# Instalar COPLA en su entorno de conda
COPY environment.yml .
RUN conda env create -f environment.yml -n copla

# Instalar MacSyFinder en un entorno separado
RUN conda env create -f macsyfinder.environment.yml -n macsyfinder

# Ejecutar conda init para configurar shells dentro del contenedor
RUN /bin/bash -c "conda init bash && \
    echo 'source /opt/conda/etc/profile.d/conda.sh' >> ~/.bashrc"

# Activar MacSyFinder y hacer instalación manual
RUN /bin/bash -c "source /opt/conda/etc/profile.d/conda.sh && \
    conda activate macsyfinder && \
    cd \${CONDA_PREFIX} && \
    wget https://github.com/gem-pasteur/macsyfinder/archive/macsyfinder-1.0.5.tar.gz && \
    tar zxf macsyfinder-1.0.5.tar.gz && \
    wget https://github.com/gem-pasteur/macsyview/archive/macsyview-1.0.1.tar.gz && \
    tar zxf macsyview-1.0.1.tar.gz && \
    mv macsyfinder-macsyfinder-1.0.5 src && \
    mv macsyview-macsyview-1.0.1 src/macsyview && \
    cd src && \
    python setup.py build && \
    python setup.py test -vv && \
    python setup.py install --prefix=\${CONDA_PREFIX} && \
    conda deactivate"

# Descargar el script ANI y moverlo a la carpeta bin
RUN wget https://github.com/lmrodriguezr/enveomics/raw/master/Scripts/ani.rb && \
    chmod +x ani.rb && \
    mv ani.rb /opt/copla/bin/

# Configurar conda para activarse automáticamente con el entorno COPLA
RUN echo "source /opt/conda/etc/profile.d/conda.sh && conda activate copla" >> ~/.bashrc

# Descargar bases de datos de COPLA
RUN bin/download_Copla_databases.sh

# Copiar el script de activación de MacSyFinder dentro del contenedor
RUN echo 'source /opt/conda/etc/profile.d/conda.sh && conda activate macsyfinder && exec "$@"' > /usr/local/bin/run_macsyfinder.sh
RUN chmod +x /usr/local/bin/run_macsyfinder.sh

# Configurar la variable de entorno para asegurar que conda puede activar entornos
ENV PATH="/opt/conda/envs/copla/bin:$PATH"

# Establecer el entorno activado por defecto
SHELL ["conda", "run", "-n", "copla", "/bin/bash", "-c"]
