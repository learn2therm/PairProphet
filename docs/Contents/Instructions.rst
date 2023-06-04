Instructions and Use Case
=============
  
| **Training**
**************** 

| PairProphet is trained on UW's Learn2Therm database (n=___ million protein pairs).

| **Use Case**
**************** 

| Assess protein pair functionality with model developed during training with **user_input.py**.

| Input: CSV file with protein pairs.
| Output: 
| CSV with:
    - protein pairs 
    - hmmer functionality (Boolean (1=True, 0=False)) 
    - structure functionality (Boolean (1=True, 0=False)) 
| .txt file with:
    - accuracy
    - mean precision
    - mean F1 score
    - mean recall

