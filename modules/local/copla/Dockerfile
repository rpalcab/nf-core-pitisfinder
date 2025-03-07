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

# Establecer el directorio de trabajo inicial en la raíz
WORKDIR /data/app/

# Clonar COPLA en un directorio dedicado (/copla)
# RUN git clone -b dev --single-branch https://github.com/rpalcab/COPLA.git /copla
COPY copla .

# Añadir copla al $PATH del sistema
RUN echo 'export PATH="$PATH:/data/app/bin"' >> ~/.bashrc

# Instalar entornos de COPLA y MacSyFinder
RUN conda env create -f /data/app/copla.environment.yml -n copla && conda env create -f /data/app/macsyfinder.environment.yml -n macsyfinder 

# Ejecutar conda init para bash y configurar el sourcing de conda.sh
RUN /bin/bash -c "conda init bash && echo 'source /opt/conda/etc/profile.d/conda.sh && conda activate copla' >> ~/.bashrc"

# Activar macsyfinder y realizar la instalación manual de MacSyFinder
RUN /bin/bash -c "source /opt/conda/etc/profile.d/conda.sh && \
    conda activate macsyfinder && \
    cd ${CONDA_PREFIX} && \
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

# Descargar el script ani.rb y moverlo a la carpeta bin de COPLA
RUN wget https://github.com/lmrodriguezr/enveomics/raw/master/Scripts/ani.rb -P /data/app/bin/ && \
    chmod +x /data/app/bin/ani.rb

# Crear un enlace simbólico para que el ejecutable se encuentre en PATH
RUN ln -s /data/app/bin/copla.py /usr/local/bin/copla && \
    ln -s /data/app/bin/download_Copla_databases.sh /usr/local/bin/download_Copla_databases && \
    ln -s /data/app/bin/generate_test_output.sh /usr/local/bin/generate_test_output && \
    ln -s /data/app/bin/post_install_test.sh /usr/local/bin/post_install_test

# Copiar el script de activación de MacSyFinder dentro del contenedor
RUN echo 'source /opt/conda/etc/profile.d/conda.sh && conda activate macsyfinder && exec "$@"' > /usr/local/bin/run_macsyfinder.sh
RUN chmod +x /usr/local/bin/run_macsyfinder.sh

# Ajustar la variable PATH para que se incluya el entorno COPLA
ENV PATH="/opt/conda/envs/copla/bin:$PATH"

# Establecer el entorno activado por defecto
SHELL ["conda", "run", "-n", "copla", "/bin/bash", "-c"]
