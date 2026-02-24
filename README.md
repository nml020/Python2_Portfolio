# NatalieLeBlanc_Python2_Portfolio
This is the portfolio of the Python 2 codes that I learned through Winter 2025-2026 in the Python 2 class (BISC 4503).








## Feature Detection
In this analysis, we took images of sunflowers for finding smaller train images in larger test images using using matching methods with heat map visualization and rectangles for specific location detection.

```python
import cv2
```


```python
import numpy as np
```


```python
import matplotlib.pyplot as plt
```


```python
%matplotlib inline
```


```python
full = cv2.imread('training_sunflower.jpg')
```


```python
full = cv2.cvtColor(full, cv2.COLOR_BGR2RGB)
```


```python
plt.imshow(full)
```




    <matplotlib.image.AxesImage at 0x7ff6c62b0e10>

<img width="441" height="474" alt="image" src="https://github.com/user-attachments/assets/580f3481-9e16-4bc8-9e16-6f58a169f79f" />



```python
test = cv2.imread('testing_sunflower.jpg')
```


```python
test = cv2.cvtColor(test, cv2.COLOR_BGR2RGB)
```


```python
plt.imshow(test)
```




    <matplotlib.image.AxesImage at 0x7ff6c621dc10>
<img width="466" height="351" alt="image" src="https://github.com/user-attachments/assets/4c6ddeb7-8e71-4159-88a5-9d4211d713c8" />




```python
print('Test image shape:', full.shape)
print('Training image shape:', test.shape)
```

    Test image shape: (1555, 1404, 3)
    Training image shape: (310, 510, 3)



```python
methods = ['cv2.TM_CCOEFF', 'cv2.TM_CCOEFF_NORMED', 'cv2.TM_CCORR', 'cv2.TM_CCORR_NORMED', 'cv2.TM_SQDIFF', 'cv2.TM_SQDIFF_NORMED']

```


```python
for m in methods:
    plt.figure()
    
    test_copy = test.copy()
    method = eval(m)
    
    res= cv2.matchTemplate(test_copy, full, method)
    
    min_val, max_val, min_loc, max_loc = cv2.minMaxLoc(res)
    
    if method in [cv2.TM_SQDIFF, cv2.TM_SQDIFF_NORMED]:
        top_left = min_loc
    else:
        top_left = max_loc
    height, width, channels = full.shape
    bottom_right = (top_left[0] + width, top_left[1] + height)
    
    cv2.rectangle(test_copy, top_left, bottom_right, (255,0,0), 10)
    
    plt.subplot(121)
    plt.imshow(res)
    plt.title("Heatmap of template matching")
    plt.subplot(122)
    plt.imshow(test_copy)
    plt.title("Detection of template")
    
    plt.suptitle(m)
    
    plt.show
    print('\n')
    print('\n')
```
<img width="659" height="895" alt="image" src="https://github.com/user-attachments/assets/ac9f6ba4-a502-4045-866d-8476ae7d6c5e" />

<img width="662" height="893" alt="image" src="https://github.com/user-attachments/assets/019cdd06-bbf9-49e8-a24f-5ed5fd717b23" />

<img width="635" height="932" alt="image" src="https://github.com/user-attachments/assets/e0d2e34b-6392-4602-a9ee-d80e3359b8f7" />

