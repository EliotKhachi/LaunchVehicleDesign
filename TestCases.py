import pandas as pd
from pptx import Presentation
from pptx.enum.shapes import MSO_SHAPE
from pptx.enum.dml import MSO_THEME_COLOR
from pptx.util import Inches
#case = {'Case': ['Step 2 %dv', 'required dv2', 'stage2\'s loss', 'staging speed', 'mu_2', 'm_p2', 'm_s2', 'm_02']}
#payload_items = ['PLF', 'Payload', 'PAF']
#df = pd.DataFrame(columns=['Item', 'Height (m)', 'Mass (kg)', 'Distance (m)', 'Moment (kg*m)', 'Thickness (m)', 'Distance from CM (m)', 'J0 (kg m^2)', 'm*CM^2 (kg m^2)', 'Jpitch/yaw', 'Jroll'], index = range(len(payload_items)))

#new_row = {'Item': 'Forward Skirt', 'Height (m)': 0, 'Mass (kg)': 0, 'Distance (m)': 0, 'Moment (kg*m)': 0, 'Thickness (m)': 0, 'Distance from CM (m)': 0, 'J0 (kg m^2)': 0, 'm*CM^2 (kg m^2)': 0, 'Jpitch/yaw': 0, 'Jroll': 0}
#df.append(new_row, ignore_index = True)

#print(df)

dict = {'Name': ['Sally', 'Jerry'],
      'Age': [25, 27],
      'Occupation': ['Therapist', 'Engineer']}
df2 = pd.DataFrame(dict)
new_row ={'Name': 'Bob'}
df2 = df2.append(new_row, ignore_index=True)

print(len(df2))
print(df2)
#df2 = df2.reindex([2,1,0])
#df2 = df2.reindex(list(range(len(df2)))[::-1])
print(df2.iloc[0:2,:])
df3 = df2.iloc[0:2,:]
print(df3)
df4 = pd.DataFrame(
      {
            'Name': NAN,
            'Age': NAN,
            'Occupation': NAN
      }
)
# print(list(range(10))[::-1])

# print(list(range(3,8)))