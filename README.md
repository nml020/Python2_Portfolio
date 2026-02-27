# NatalieLeBlanc_Python2_Portfolio
This is the portfolio of the Python 2 codes that I learned through Winter 2025-2026 in the Python 2 class (BISC 4503).


## Sequence Objects


## Sequence Annotations


## Sequence I/O


## Multiple Sequence Alignment 



## Blast 


## Challenge 1


## Open CV
### Open CV 1: "OpenCVBasics"
In this analysis, we used an image, particularly of a sloth, to learn the basics of image processing with modifications including color conservion and image rotation. 

```python
import numpy as np
import matplotlib.pyplot as plt
%matplotlib inline
```


```python
import cv2
```


```python
img = cv2.imread("sloth.jpeg")
```


```python
type(img)
```




    numpy.ndarray




```python
img_wrong = cv2.imread('wrong/path/doesnt/abcdegh.jpg')
```


```python
type(img_wrong)
```




    NoneType




```python
plt.imshow(img)
```




    <matplotlib.image.AxesImage at 0x7efcdf41fd50>

![Uploading image.png…]()




```python
fix_img = cv2.cvtColor(img, cv2.COLOR_BGR2RGB)
```


```python
plt.imshow(fix_img)
```




    <matplotlib.image.AxesImage at 0x7efcdf3d6dd0>


![Uploading image.png…]()



```python
img_gray = cv2.imread("sloth.jpeg", cv2.IMREAD_GRAYSCALE)
img_gray.shape
```




    (190, 266)




```python
plt.imshow(img_gray)
```




    <matplotlib.image.AxesImage at 0x7efcddb49b10>


![Uploading image.png…]()




```python
plt.imshow(img_gray, cmap = "gray")
```




    <matplotlib.image.AxesImage at 0x7efcddabc710>


![Uploading image.png…]()




```python
fix_img.shape
```




    (190, 266, 3)




```python
new_img = cv2.resize(fix_img,(1000,400))
plt.imshow(new_img)
```




    <matplotlib.image.AxesImage at 0x7efcdc224f10>

![Uploading image.png…]()




```python
new_img.shape
```




    (400, 1000, 3)




```python
w_ratio = 0.5
h_ratio = 0.5

new_img = cv2.resize(fix_img, (0,0), fix_img, w_ratio, h_ratio)
```


```python
plt.imshow(new_img)
```




    <matplotlib.image.AxesImage at 0x7efcdc207a10>


![Uploading image.png…]()




```python
new_img.shape
```




    (95, 133, 3)




```python
flip_img = cv2.flip(fix_img, 0)
plt.imshow(flip_img)
```




    <matplotlib.image.AxesImage at 0x7efcdc174690>


![Uploading image.png…]()




```python
flip_img2 = cv2.flip(fix_img, -1)
plt.imshow(flip_img2)
```




    <matplotlib.image.AxesImage at 0x7efcdc0e3790>

![Uploading image.png…]()




```python
type(fix_img)
```




    numpy.ndarray




```python
cv2.imwrite('sloth_fixed_image.jpeg', flip_img)
```




    True



## Aspect Detection
### Corner Detection
In this analysis, we used to chessboard images to detect and visualize their corner features using Harris and Shi-Tomasi corner detection methods. 
```python
import cv2
import numpy as np
import matplotlib.pyplot as plt
%matplotlib inline
```


```python
flat_chess = cv2.imread('chessboard_green.png')
flat_chess = cv2.cvtColor(flat_chess, cv2.COLOR_BGR2RGB)
plt.imshow(flat_chess)
```




    <matplotlib.image.AxesImage at 0x7f071da3cf90>

<img width="450" height="403" alt="image" src="https://github.com/user-attachments/assets/f64e5d4b-1975-404f-9ef1-d222fca6965c" />




```python
gray_flat_chess = cv2.cvtColor(flat_chess, cv2.COLOR_BGR2GRAY)
plt.imshow(gray_flat_chess, cmap = "gray")
```




    <matplotlib.image.AxesImage at 0x7f071d9aaa10>


<img width="457" height="411" alt="image" src="https://github.com/user-attachments/assets/e25f604c-e734-435b-bdb1-a1e7eb4cc2a2" />




```python
real_chess = cv2.imread("chessboard.jpg")
real_chess = cv2.cvtColor(real_chess, cv2.COLOR_BGR2RGB)
```


```python
plt.imshow(real_chess)
```




    <matplotlib.image.AxesImage at 0x7f071c101c90>

<img width="643" height="413" alt="image" src="https://github.com/user-attachments/assets/6e1efa85-5e6d-4ac0-9d0d-4d25a91c80d3" />




```python
gray_real_chess = cv2. cvtColor(real_chess, cv2.COLOR_BGR2GRAY)
plt.imshow(gray_real_chess, cmap = 'gray')
```




    <matplotlib.image.AxesImage at 0x7f071c0950d0>


<img width="662" height="427" alt="image" src="https://github.com/user-attachments/assets/72681734-abe9-4e3c-8110-fa6cbf4e38ab" />




```python
gray = np.float32(gray_flat_chess)
dst = cv2.cornerHarris(src = gray, blockSize = 2, ksize = 3, k = 0.04)

dst = cv2.dilate(dst, None)
```


```python
flat_chess[dst>0.01*dst.max()] = [255,0,0]

plt.imshow(flat_chess)
```




    <matplotlib.image.AxesImage at 0x7f071c0df9d0>

<img width="466" height="421" alt="image" src="https://github.com/user-attachments/assets/8babc0e1-6b37-40c9-9165-d10a25d06e19" />




```python
gray = np.float32(gray_real_chess)
dst = cv2.cornerHarris(src = gray, blockSize = 2, ksize = 3, k = 0.04)
dst = cv2.dilate(dst, None)

real_chess[dst>0.01*dst.max()] = [255, 0, 0]

plt.imshow(real_chess)
```




    <matplotlib.image.AxesImage at 0x7f0715f3d650>


<img width="642" height="407" alt="image" src="https://github.com/user-attachments/assets/1cdb3322-131c-4b88-9237-cad8a9621a06" />




```python
#Shi-Tomasi Corner Detection

corners = cv2.goodFeaturesToTrack(gray_flat_chess, 64, 0.01, 10)
```


```python
corners = np.int0(corners)

for i in corners:
    x,y = i.ravel()
    cv2.circle(flat_chess, (x,y,), 3, (255,0,0), -1)
    
plt.imshow(flat_chess)
```




    <matplotlib.image.AxesImage at 0x7f0715ec0b90>


<img width="446" height="405" alt="image" src="https://github.com/user-attachments/assets/2bbc80e4-33ed-4fc6-8e2b-6e3a3c0f616e" />



```python
corners = cv2.goodFeaturesToTrack(gray_real_chess, 100, 0.01, 10)

corners = np.int0(corners)

for i in corners:
    x,y = i.ravel()
    cv2.circle(real_chess, (x,y), 3, (0,255,0), -1)
    
plt.imshow(real_chess)
```




    <matplotlib.image.AxesImage at 0x7f0715e97b50>

<img width="712" height="407" alt="image" src="https://github.com/user-attachments/assets/7ae4a27b-be2b-43d7-9ccd-c35ca37a7466" />




### Edge Detection


## Feature Detection
### Feature Matches
In this analysis, we took images of Apple Jacks and a variety of cereal and used ORB and SIFT to detect and match, using BF and FLANN matching, within one another with visualizations through line connections.

```python
import cv2
import numpy as np
import matplotlib.pyplot as plt
%matplotlib inline
```


```python
def display(img,cmap = 'gray'):
    fig = plt.figure(figsize =(12, 10))
    ax = fig.add_subplot(111)
    ax.imshow(img,cmap = 'gray')
```


```python
apple_jacks = cv2.imread("apple_jacks.jpg", 0)
display(apple_jacks)
```
<img width="609" height="930" alt="image" src="https://github.com/user-attachments/assets/dd5f541a-3d30-4f4d-8249-94a392a86567" />




```python
cereals = cv2.imread('all_cereal.jpg', 0)
display(cereals)
```
<img width="1168" height="857" alt="image" src="https://github.com/user-attachments/assets/28e6a75d-d889-4ac8-8349-e431ce530640" />




```python
orb = cv2.ORB_create()

kp1,des1 = orb.detectAndCompute(apple_jacks, mask=None)
kp2,des2 = orb.detectAndCompute(cereals, mask=None)
```


```python
bf = cv2.BFMatcher(cv2.NORM_HAMMING, crossCheck = True)
matches = bf.match(des1, des2)
```


```python
matches = sorted(matches, key = lambda x:x.distance)
```


```python
apple_jacks_matches = cv2.drawMatches(apple_jacks, kp1, cereals, kp2, matches[:25], None, flags = 2)
```


```python
display(apple_jacks_matches)
```

<img width="1160" height="714" alt="image" src="https://github.com/user-attachments/assets/2404ff72-3e52-41d4-8eec-b7001ed0cce8" />



```python
sift = cv2.SIFT_create()
```


```python
kp1 = sift.detect(apple_jacks, None)
kp1, des1 = sift.compute(apple_jacks, kp1)

kp2 = sift.detect(cereals, None)
kp2, des2 = sift.compute(cereals, kp2)
```


```python
bf = cv2.BFMatcher()
macthes = bf.knnMatch(des1, des2, k=2)
```


```python
good = []

for match1, match2 in macthes:
    if match1.distance < 0.75*match2.distance:
        good.append([match1])
```


```python
print('Length of total matches:', len(macthes))
print('Length of good matches:', len(good))
```

    Length of total matches: 4316
    Length of good matches: 120



```python
sift_matches = cv2.drawMatchesKnn(apple_jacks, kp1, cereals, kp2, good, None, flags = 2)
display(sift_matches)
```

<img width="1157" height="708" alt="image" src="https://github.com/user-attachments/assets/e2b929ff-a0ba-47ef-a11c-539350c35a81" />



```python
sift = cv2.SIFT_create()

kp1 = sift.detect(apple_jacks, None)
kp1, des1 = sift.compute(apple_jacks, kp1)

kp2 = sift.detect(cereals, None)
kp2, des2 = sift.compute(cereals, kp2)
```


```python
flann_index_KDtree = 0
index_params = dict(algorithm=flann_index_KDtree, trees = 5)
search_params = dict(checks=50)
```


```python
flann = cv2.FlannBasedMatcher(index_params, search_params)

matches = flann.knnMatch(des1, des2, k=2)

good = []

for match1, match2, in matches:
    if match1.distance < 0.75*match2.distance:
        good.append([match1])
```


```python
flann_matches = cv2.drawMatchesKnn(apple_jacks, kp1, cereals, kp2, good, None, flags = 0)
display(flann_matches)
```

<img width="1271" height="721" alt="image" src="https://github.com/user-attachments/assets/6263de34-2246-4797-b8c0-96654915177e" />



```python
sift = cv2.SIFT_create()

kp1 = sift.detect(apple_jacks, None)
kp1, des1 = sift.compute(apple_jacks, kp1)

kp2 = sift.detect(cereals, None)
kp2, des2 = sift.compute(cereals, kp2)
```


```python
flann_index_KDtree = 0
index_params = dict(algorithm = flann_index_KDtree, trees = 5)
search_param = dict(checks = 50)
```


```python
flann = cv2.FlannBasedMatcher(index_params, search_params)

matches = flann.knnMatch(des1, des2, k = 2)
```


```python
matchesMask = [[0,0] for i in range(len(matches))]
```


```python
for i, (match1, match2) in enumerate(matches):
    if match1.distance <0.75*match2.distance:
        matchesMask[i] = [1,0]

draw_params = dict(matchColor = (0,255,0),
                  singlePointColor = (255,0,0),
                  matchesMask = matchesMask,
                  flags = 0)
```


```python
flann_matches = cv2.drawMatchesKnn(apple_jacks, kp1, cereals, kp2, matches, None, **draw_params)

display(flann_matches)
```
<img width="1152" height="716" alt="image" src="https://github.com/user-attachments/assets/7d4fe00e-5178-424e-99e2-5f41adc29208" />



### Object Detection
In this analysis, we took images of sunflowers for finding smaller train images in larger test images using using matching methods with heat map visualization and rectangles for specific location detection.

```python
# Import libraries
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
# Load training image 
full = cv2.imread('training_sunflower.jpg')
```


```python
# Convert image to color
full = cv2.cvtColor(full, cv2.COLOR_BGR2RGB)
```


```python
# Display image
plt.imshow(full)
```




    <matplotlib.image.AxesImage at 0x7ff6c62b0e10>

<img width="441" height="474" alt="image" src="https://github.com/user-attachments/assets/580f3481-9e16-4bc8-9e16-6f58a169f79f" />



```python
# Load testing image
test = cv2.imread('testing_sunflower.jpg')
```


```python
# Convert image to color
test = cv2.cvtColor(test, cv2.COLOR_BGR2RGB)
```


```python
# Display image
plt.imshow(test)
```




    <matplotlib.image.AxesImage at 0x7ff6c621dc10>
<img width="466" height="351" alt="image" src="https://github.com/user-attachments/assets/4c6ddeb7-8e71-4159-88a5-9d4211d713c8" />




```python
# Print image dimensions
print('Test image shape:', full.shape)
print('Training image shape:', test.shape)
```

    Test image shape: (1555, 1404, 3)
    Training image shape: (310, 510, 3)



```python
# Define matching methods
methods = ['cv2.TM_CCOEFF', 'cv2.TM_CCOEFF_NORMED', 'cv2.TM_CCORR', 'cv2.TM_CCORR_NORMED', 'cv2.TM_SQDIFF', 'cv2.TM_SQDIFF_NORMED']

```


```python
# Instruction loop for each macthing method: copy, method, template, location, box detection
# Heatmap display with detection results
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

