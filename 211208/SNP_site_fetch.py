from selenium import webdriver
from selenium.webdriver import FirefoxOptions
from selenium.webdriver.common.by import By
from time import sleep

fo = open("geneweblist.tsv", "r")

for line in fo:
    profile = webdriver.FirefoxProfile()
    profile.set_preference('browser.download.folderList', 2)
    profile.set_preference('browser.download.manager.showWhenStarting', False)
    profile.set_preference('browser.download.dir', '~/temp/')
    profile.set_preference('browser.helperApps.neverAsk.saveToDisk', 'text/csv')
    opts = FirefoxOptions()
    opts.add_argument("--headless")
    browser = webdriver.Firefox(firefox_profile=profile, options=opts)
    browser.get(line)
    sleep(4.5)
    print(browser.title)
    button = browser.find_elements(By.CLASS_NAME, 'Button__BaseButton-sc-1eobygi-0.Button-sc-1eobygi-1.indcWT')
    for i in button:
        if i.text == 'Export variants to CSV':
            i.click()
    sleep(1.5)
    browser.close()

fo.close()
