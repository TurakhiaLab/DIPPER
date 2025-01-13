import re

def remove_braces_content(text):
    return re.sub(r'\{.*?\}', '', text)

input_text = input()

output_text = remove_braces_content(input_text)

print(output_text)