# liteRF
Световые поля

# Устанока
1. Клонировать репозиторий:X
   * git clone https://github.com/Lorent1/liteRF
   * cd liteRF 
2. Сабмодули:
   * git submodule init && git submodule update 
3. Установить модель (модель 2 - 256x256x256):
   * https://drive.google.com/file/d/1hhns6AGGUv28U9DQVfnYYTRjTdp60NHT/view?usp=drive_link

# Сборка (CPU):
Сборка, используя Cmake
   * mkdir build-cpu && cd cmake-cpu
   * cmake -DCMAKE_BUILD_TYPE=Release ..
   * make -j 8

# Сборка (GPU) (ВНИМАНИЕ! Версия GPU не будет работать корректно, проблемы с копированием сетки на GPU, не хватило времени): 
Сборка, используя Cmake с флагом 'USE_VULKAN' == 'ON':
   * mkdir build-gpu && cd build-gpu
   * cmake -DCMAKE_BUILD_TYPE=Release -DUSE_VULKAN=ON ..
   * make -j 8

# Время работы

Результатом программы является файлы вида `out_cpu_{i}.bmp`

Все расчеты проводились на ноутбуке - HP Victus / `Ryzen 5 5600H` / `GeForce RTX 3050Ti 4GB` / 16GB

За основу была взята модель размеров `256^3`, в таблице указано время на отрисовку одного кадра в разрешении `512x512`

| Вычислитель/Разрешение      | 512x512      |
| --------------------------- | -------------|
| CPU (1 поток)               | 25.722 c     |
| CPU (многопоточный)         | 4.227 с      | 
| GPU (без учета копирования) | -            |
| GPU (с учетом копирования)  | -            |
| GPU (время копирования)     | -            | 
